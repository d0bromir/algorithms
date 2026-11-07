#!/usr/bin/env python3
"""Static validation of the WES workflows against the real Galaxy tool wrappers (no Galaxy server needed).

For every tool step (including the embedded subworkflow) it checks that
  * the Tool Shed tool id/version exists in the wrapper sources (galaxyproject/tools-iuc),
  * tool_state + input_connections are a valid "workflow_step_linked" state according to Galaxy's own
    parameter models (galaxy-tool-util): unknown keys, invalid select/conditional values, wrong types,
  * every connection targets an existing data/collection (or, for workflow parameters, scalar) parameter
    of the selected conditional case, and comes from an existing output of the source step,
  * the datatype of every connected output is accepted by the receiving parameter (Galaxy datatype
    registry from galaxy-data, including implicit converters).

    pip install -r requirements.txt
    git clone --depth 1 https://github.com/galaxyproject/tools-iuc.git /tmp/tools-iuc
    python validate_workflows.py --tools-iuc /tmp/tools-iuc            # validate ../*.ga
    python validate_workflows.py --tools-iuc /tmp/tools-iuc --self-test  # also prove that bad states are caught
"""
import argparse
import copy
import glob
import json
import os
import re
import sys

from galaxy.tool_util.parameters import input_models_for_tool_source, validate_workflow_step_linked
from galaxy.tool_util.parser import get_tool_source

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL_SHED = "toolshed.g2.bx.psu.edu"

# Tool Shed repository (owner, name) -> wrapper directory in tools-iuc and the tool XML files used.
REPOSITORIES = [
    ("devteam", "fastqc", "tools/fastqc", ["rgFastQC.xml"]),
    ("iuc", "fastp", "tools/fastp", ["fastp.xml"]),
    ("devteam", "bwa", "tools/bwa", ["bwa-mem.xml"]),
    ("devteam", "picard", "tools/picard",
     ["picard_MarkDuplicates.xml", "picard_BedToIntervalList.xml", "picard_CollectHsMetrics.xml"]),
    ("devteam", "samtools_stats", "tool_collections/samtools/samtools_stats", ["samtools_stats.xml"]),
    ("iuc", "samtools_faidx", "tool_collections/samtools/samtools_faidx", ["samtools_faidx.xml"]),
    ("iuc", "mosdepth", "tools/mosdepth", ["mosdepth.xml"]),
    ("devteam", "freebayes", "tools/freebayes", ["freebayes.xml"]),
    ("iuc", "bedtools", "tools/bedtools", ["slopBed.xml"]),
    ("iuc", "bcftools_norm", "tools/bcftools", ["bcftools_norm.xml"]),
    ("iuc", "bcftools_filter", "tools/bcftools", ["bcftools_filter.xml"]),
    ("iuc", "bcftools_annotate", "tools/bcftools", ["bcftools_annotate.xml"]),
    ("iuc", "bcftools_stats", "tools/bcftools", ["bcftools_stats.xml"]),
    ("iuc", "snpeff", "tool_collections/snpeff", ["snpEff.xml"]),
    ("iuc", "snpsift", "tool_collections/snpsift/snpsift", ["snpSift_filter.xml", "snpSift_extractFields.xml"]),
    ("iuc", "multiqc", "tools/multiqc", ["multiqc.xml"]),
]

# Galaxy built-in tools: flat inputs (type, formats, multiple) and outputs.
BUILTIN_TOOLS = {
    "__FLATTEN__": {
        "inputs": {"input": ("data_collection", None, False), "join_identifier": ("select", None, False)},
        "outputs": {"output": None},
    },
}


class Validator:
    def __init__(self, tool_index, registry):
        self.tool_index = tool_index
        self.registry = registry
        self.errors = []
        self._sources = {}

    def err(self, msg):
        self.errors.append(msg)
        print("  ERROR:", msg)

    def tool_source(self, tid):
        if tid not in self._sources:
            self._sources[tid] = get_tool_source(self.tool_index[tid])
        return self._sources[tid]

    # ---------------------------------------------------------------------------------------- outputs
    def step_outputs(self, step):
        """{output_name: datatype or None} for any step type (post-job datatype changes applied)."""
        t = step["type"]
        if t in ("data_input", "data_collection_input"):
            fmts = json.loads(step.get("tool_state") or "{}").get("format") or []
            return {"output": fmts[0] if len(fmts) == 1 else None}
        if t == "parameter_input":
            return {"output": None}
        if t == "subworkflow":
            result = {}
            for s in step["subworkflow"]["steps"].values():
                for wo in s.get("workflow_outputs") or []:
                    result[wo["label"]] = self.step_outputs(s).get(wo["output_name"])
            return result
        tid = step["tool_id"]
        if tid in BUILTIN_TOOLS:
            return dict(BUILTIN_TOOLS[tid]["outputs"])
        outs, colls = self.tool_source(tid).parse_outputs(None)
        result = {name: getattr(o, "format", None) for name, o in outs.items()}
        result.update({name: None for name in colls})
        for pja in (step.get("post_job_actions") or {}).values():
            if pja["action_type"] == "ChangeDatatypeAction":
                result[pja["output_name"]] = pja["action_arguments"]["newtype"]
        return result

    # -------------------------------------------------------------------------------------- datatypes
    def accepts(self, src_fmt, accepted):
        reg = self.registry
        if reg is None or src_fmt is None or not accepted:
            return True
        src = reg.get_datatype_by_extension(src_fmt)
        if src is None:
            return True
        for a in (x.strip() for x in accepted):
            if a in ("data", "input"):
                return True
            target = reg.get_datatype_by_extension(a)
            if target is not None and isinstance(src, type(target)):
                return True
            if a in reg.datatype_converters.get(src_fmt, {}):
                return True
        return False

    # ------------------------------------------------------------------------------------ validation
    def validate(self, wf, indent=""):
        steps = wf["steps"]
        print(f"{indent}Workflow: {wf.get('name')} ({len(steps)} steps)")
        for sid, step in sorted(steps.items(), key=lambda x: int(x[0])):
            label = step.get("label") or sid
            where = f"{wf.get('name')} / step {sid} ({label})"
            if step["type"] == "subworkflow":
                print(f"{indent}- subworkflow {label}")
                self.validate_subworkflow_step(step, steps, where)
                self.validate(step["subworkflow"], indent + "    ")
            elif step["type"] == "tool":
                print(f"{indent}- {label}: {step['tool_id']}")
                self.validate_tool_step(step, steps, where)

    def validate_subworkflow_step(self, step, steps, where):
        if "subworkflow" not in step:
            self.err(f"{where}: subworkflow is not embedded")
            return
        sub_inputs = {s["label"]: s for s in step["subworkflow"]["steps"].values()
                      if s["type"] in ("data_input", "data_collection_input", "parameter_input")}
        connections = step.get("input_connections") or {}
        for name, conns in connections.items():
            target = sub_inputs.get(name)
            if target is None:
                self.err(f"{where}: connection to unknown subworkflow input {name!r}")
                continue
            for c in conns:
                outs = self.step_outputs(steps[str(c["id"])])
                if c["output_name"] not in outs:
                    self.err(f"{where}: source step {c['id']} has no output {c['output_name']!r}")
                elif target["type"] == "data_input":
                    accepted = json.loads(target["tool_state"]).get("format")
                    if not self.accepts(outs[c["output_name"]], accepted):
                        self.err(f"{where}: {outs[c['output_name']]} not accepted by subworkflow input {name}")
        for name, s in sub_inputs.items():
            st = json.loads(s["tool_state"])
            if name not in connections and not st.get("optional") and st.get("default") is None:
                self.err(f"{where}: required subworkflow input {name!r} is not connected")

    def validate_tool_step(self, step, steps, where):
        tid = step["tool_id"]
        if tid in BUILTIN_TOOLS:
            flat = {k: [((),) + v] for k, v in BUILTIN_TOOLS[tid]["inputs"].items()}
            bundle = None
        elif tid not in self.tool_index:
            self.err(f"{where}: tool {tid} not found in the wrapper sources")
            return
        else:
            ts = self.tool_source(tid)
            if ts.parse_version() != step.get("tool_version"):
                self.err(f"{where}: tool_version {step.get('tool_version')} != wrapper {ts.parse_version()}")
            bundle = input_models_for_tool_source(ts)
            flat = {}
            for page in ts.parse_input_pages().page_sources:
                flatten(page.parse_input_sources(), "", flat, ())
        state = strip_private(json.loads(step["tool_state"]))
        for name, conns in (step.get("input_connections") or {}).items():
            if name == "when":
                continue
            resolved = resolve_param(flat, name, state)
            if resolved is None:
                self.err(f"{where}: connection to unknown or inactive parameter {name!r}")
                continue
            ptype, fmts, _multiple = resolved
            for c in conns:
                src = steps[str(c["id"])]
                outs = self.step_outputs(src)
                if c["output_name"] not in outs:
                    self.err(f"{where}: source step {c['id']} ({src.get('label')}) has no output {c['output_name']!r}")
                    continue
                if src["type"] == "parameter_input":
                    ptype_src = json.loads(src["tool_state"]).get("parameter_type")
                    compatible = {"text": {"text", "select"}, "integer": {"integer", "float", "text"},
                                  "float": {"float", "text"}, "boolean": {"boolean"}}
                    if ptype not in compatible.get(ptype_src, set()):
                        self.err(f"{where}: {ptype_src} workflow parameter connected to {ptype} parameter {name}")
                else:
                    if ptype not in ("data", "data_collection"):
                        self.err(f"{where}: dataset connected to non-data parameter {name} ({ptype})")
                    elif ptype == "data" and not self.accepts(outs[c["output_name"]], fmts):
                        self.err(f"{where}: {outs[c['output_name']]} from {src.get('label')}/{c['output_name']} "
                                 f"not accepted by {name} {fmts}")
            set_connected(state, name)
        if bundle is not None:
            try:
                validate_workflow_step_linked(bundle, state)
            except Exception as e:  # pydantic ValidationError / RequestParameterInvalidException
                self.err(f"{where}: invalid tool state: {e}")
        outs = self.step_outputs({**step, "post_job_actions": {}})
        for pja in (step.get("post_job_actions") or {}).values():
            if pja.get("output_name") and pja["output_name"] not in outs:
                self.err(f"{where}: post job action on unknown output {pja['output_name']}")
        for wo in step.get("workflow_outputs") or []:
            if wo["output_name"] not in outs:
                self.err(f"{where}: workflow output on unknown tool output {wo['output_name']}")


# ------------------------------------------------------------------------------------------- helpers
def flatten(inputs, prefix, out, conds):
    """path -> [(conditions, type, formats, multiple)], conditions = ((test_path, case_value), ...)."""
    for inp in inputs:
        t = inp.parse_input_type()
        name = inp.parse_name()
        if t == "param":
            fmts = inp.get("format")
            out.setdefault(prefix + name, []).append(
                (conds, inp.get("type"), fmts.split(",") if fmts else None, inp.get_bool("multiple", False)))
        elif t == "conditional":
            test = inp.parse_test_input_source()
            test_path = f"{prefix}{name}|{test.parse_name()}"
            out.setdefault(test_path, []).append((conds, test.get("type"), None, False))
            for value, case in inp.parse_when_input_sources():
                flatten(case.parse_input_sources(), f"{prefix}{name}|", out, conds + ((test_path, value),))
        elif t == "section":
            flatten(inp.parse_nested_inputs_source().parse_input_sources(), f"{prefix}{name}|", out, conds)
        elif t == "repeat":
            for i in range(10):
                flatten(inp.parse_nested_inputs_source().parse_input_sources(), f"{prefix}{name}_{i}|", out, conds)
    return out


def state_value(state, path):
    node = state
    for part in path.split("|"):
        m = re.match(r"^(.*)_(\d+)$", part)
        if m and isinstance(node, dict) and isinstance(node.get(m.group(1)), list):
            lst, idx = node[m.group(1)], int(m.group(2))
            if idx >= len(lst):
                return None
            node = lst[idx]
        elif isinstance(node, dict) and part in node:
            node = node[part]
        else:
            return None
    return node


def resolve_param(flat, path, state):
    candidates = flat.get(path, [])
    for conds, *definition in candidates:
        if all(str(state_value(state, tp)) == str(v) for tp, v in conds):
            return tuple(definition)
    return tuple(candidates[0][1:]) if len(candidates) == 1 else None


def set_connected(state, path):
    parts = path.split("|")
    node = state
    for i, part in enumerate(parts):
        last = i == len(parts) - 1
        m = re.match(r"^(.*)_(\d+)$", part)
        if m and isinstance(node.get(m.group(1)), list):
            lst, idx = node[m.group(1)], int(m.group(2))
            while len(lst) <= idx:
                lst.append({})
            if last:
                lst[idx] = {"__class__": "ConnectedValue"}
            else:
                node = lst[idx]
        elif last:
            node[part] = {"__class__": "ConnectedValue"}
        else:
            if not isinstance(node.get(part), dict):
                node[part] = {}
            node = node[part]


def strip_private(d):
    if isinstance(d, dict):
        return {k: strip_private(v) for k, v in d.items() if not k.startswith("__") or k == "__class__"}
    if isinstance(d, list):
        return [strip_private(v) for v in d]
    return d


def build_tool_index(tools_iuc):
    index = {}
    for owner, repo, directory, xmls in REPOSITORIES:
        for xml in xmls:
            path = os.path.join(tools_iuc, directory, xml)
            ts = get_tool_source(path)
            index[f"{TOOL_SHED}/repos/{owner}/{repo}/{ts.parse_id()}/{ts.parse_version()}"] = path
    return index


def load_registry():
    import galaxy.tool_util.linters as linters
    from galaxy.datatypes.registry import Registry

    registry = Registry()
    registry.load_datatypes(root_dir=".", config=os.path.join(os.path.dirname(linters.__file__),
                                                              "datatypes_conf.xml.sample"))
    return registry


def self_test(tool_index, registry, vcf_workflow):
    """Inject one fault per check into the VCF workflow and require every one to be reported."""
    wf = copy.deepcopy(vcf_workflow)
    by_label = {s["label"]: s for s in wf["steps"].values()}
    ids = {s["label"]: k for k, s in wf["steps"].items()}

    def edit_state(label, fn):
        state = json.loads(by_label[label]["tool_state"])
        fn(state)
        by_label[label]["tool_state"] = json.dumps(state)

    edit_state("normalize", lambda s: s.update(check_ref="zz"))
    edit_state("prioritize", lambda s: s.update(bogus_param=1))
    edit_state("target_filter", lambda s: s["sec_restrict"]["regions"].update(regions_src="nope"))
    conns = by_label["candidate_table"]["input_connections"]
    conns["inputt"] = conns.pop("input")
    by_label["annotated_table"]["input_connections"]["input"][0]["output_name"] = "no_such_output"
    by_label["variant_stats"]["input_connections"]["input_file"][0] = {
        "id": int(ids["functional_annotation"]), "output_name": "statsFile"}  # html into a VCF parameter
    expected = ["(normalize)", "(target_filter)", "(prioritize)", "inputt", "no_such_output", "format html"]
    v = Validator(tool_index, registry)
    print("\nSelf-test: validating a deliberately broken copy of the VCF workflow")
    v.validate(wf)
    missing = [e for e in expected if not any(e.replace("format ", "") in msg for msg in v.errors)]
    if missing:
        print("SELF-TEST FAILED: faults not detected:", missing)
        return False
    print(f"Self-test passed: all {len(expected)} injected faults were reported ({len(v.errors)} errors).")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tools-iuc", required=True, help="checkout of https://github.com/galaxyproject/tools-iuc")
    parser.add_argument("--self-test", action="store_true", help="also check that injected faults are detected")
    parser.add_argument("workflows", nargs="*", default=sorted(glob.glob(os.path.join(HERE, "..", "*.ga"))))
    args = parser.parse_args()

    tool_index = build_tool_index(args.tools_iuc)
    registry = load_registry()
    validator = Validator(tool_index, registry)
    for path in args.workflows:
        print("=" * 100)
        print(os.path.relpath(path))
        validator.validate(json.load(open(path)))
    print("=" * 100)
    print(f"{len(validator.errors)} error(s)")
    ok = not validator.errors
    if args.self_test:
        vcf = [p for p in args.workflows if p.endswith("wes-tertiary-from-vcf.ga")]
        ok = self_test(tool_index, registry, json.load(open(vcf[0]))) and ok
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
