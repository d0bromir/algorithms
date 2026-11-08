#!/usr/bin/env python3
"""Build the importable Galaxy workflows (.ga) from the gxformat2 sources in ./src.

The .ga files are already committed; run this only after editing src/*.gxwf.yml.
Requires Python >= 3.10 and gxformat2 >= 0.29:

    python3 -m venv .venv && .venv/bin/pip install "gxformat2>=0.29"
    .venv/bin/python build.py

* ``run: {"@import": "<file>.gxwf.yml"}`` in a step is replaced by the content of that file, so the
  FASTQ workflow carries the VCF workflow as an embedded (self-contained) subworkflow.
* UUIDs are derived from the workflow/step labels, so rebuilding an unchanged source produces
  an identical .ga file. They are hashed like uuid5 but stamped as version 4, because Galaxy's API
  validates workflow UUIDs as UUID4 and returns HTTP 500 when opening a workflow with a uuid5.
* Tool states are completed from src/galaxy-tool-states.json: the states Galaxy computes when it loads
  each step, with every parameter the sources leave out set to its default. Without them Galaxy fills
  the defaults itself and shows "may contain tools which have changed" on every import. A recorded
  state is used only if the tool version matches and it agrees with every value the source sets;
  otherwise the source state is kept and a warning asks to refresh (tests/refresh_tool_states.py).
"""
import copy
import json
import os
import sys
import uuid

try:
    from gxformat2.converter import python_to_workflow
    from gxformat2.yaml import ordered_load
except ImportError:
    sys.exit(
        "build.py needs the 'gxformat2' package (Python >= 3.10):\n"
        "    python3 -m venv .venv && .venv/bin/pip install 'gxformat2>=0.29'\n"
        "    .venv/bin/python build.py\n"
        "You do not need to build anything to use the workflows: import the committed .ga files into Galaxy."
    )

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
NAMESPACE = uuid.UUID("6b1d6c0e-3f39-4a8a-9d76-6c1f0f6f2f3a")

WORKFLOWS = {
    "wes-tertiary-from-vcf.gxwf.yml": "wes-tertiary-from-vcf.ga",
    "wes-tertiary-from-fastq.gxwf.yml": "wes-tertiary-from-fastq.ga",
}


def load(name):
    with open(os.path.join(SRC, name)) as f:
        return ordered_load(f)


def inline_imports(wf):
    steps = wf.get("steps", {})
    items = steps.items() if isinstance(steps, dict) else enumerate(steps)
    for _key, step in items:
        run = step.get("run")
        if isinstance(run, dict) and "@import" in run:
            step["run"] = inline_imports(copy.deepcopy(load(run["@import"])))
        elif isinstance(run, dict):
            inline_imports(run)
    return wf


def stable_uuid4(name):
    """Deterministic UUID with version 4 bits (Galaxy rejects other versions)."""
    return str(uuid.UUID(bytes=uuid.uuid5(NAMESPACE, name).bytes, version=4))


def pin_uuids(native, scope):
    native["uuid"] = stable_uuid4(scope)
    for key, step in native["steps"].items():
        step_scope = f"{scope}/{step.get('label') or key}"
        if "uuid" in step:
            step["uuid"] = stable_uuid4(step_scope)
        for i, wo in enumerate(step.get("workflow_outputs") or []):
            if "uuid" in wo:
                wo["uuid"] = stable_uuid4(f"{step_scope}/output/{wo.get('label') or i}")
        if isinstance(step.get("subworkflow"), dict):
            pin_uuids(step["subworkflow"], step_scope)


def load_recorded_states():
    path = os.path.join(SRC, "galaxy-tool-states.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def _scalar(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    return None if v is None else str(v)


def agrees(source, recorded, path=""):
    """Every value the source sets must equal the recorded one; returns the first mismatch or None."""
    if isinstance(source, dict):
        if not isinstance(recorded, dict):
            return f"{path}: {source!r} vs {recorded!r}"
        for k, v in source.items():
            if k not in recorded:
                return f"{path}/{k}: missing in recorded state"
            bad = agrees(v, recorded[k], f"{path}/{k}")
            if bad:
                return bad
        return None
    if isinstance(source, list) or isinstance(recorded, list):  # repeats: same length, item by item
        if not (isinstance(source, list) and isinstance(recorded, list)) or len(source) != len(recorded):
            return f"{path}: {source!r} vs {recorded!r}"
        for i, (s, r) in enumerate(zip(source, recorded)):
            bad = agrees(s, r, f"{path}[{i}]")
            if bad:
                return bad
        return None
    return None if _scalar(source) == _scalar(recorded) else f"{path}: {source!r} vs {recorded!r}"


def galaxy_encode(state):
    """Encode like Galaxy's own export: top-level values JSON strings, __page__ etc. left as they are."""
    return json.dumps({k: (v if k.startswith("__") else json.dumps(v)) for k, v in state.items()})


def complete_tool_states(native, scope, recorded, report):
    for key, step in native["steps"].items():
        step_scope = f"{scope}/{step.get('label') or key}"
        if step.get("type") == "tool":
            rec = recorded.get(step_scope)
            if rec is None:
                report.append(f"no recorded state for {step_scope}")
                continue
            if (rec["tool_id"], rec["tool_version"]) != (step["tool_id"], step["tool_version"]):
                report.append(f"{step_scope}: recorded for {rec['tool_version']}, step uses {step['tool_version']}")
                continue
            source = json.loads(step["tool_state"])
            bad = agrees({k: v for k, v in source.items() if not k.startswith("__")}, rec["tool_state"])
            if bad:
                report.append(f"{step_scope}: source differs from recorded state at {bad}")
                continue
            step["tool_state"] = galaxy_encode(rec["tool_state"])
        if isinstance(step.get("subworkflow"), dict):
            complete_tool_states(step["subworkflow"], step_scope, recorded, report)


def main():
    recorded = load_recorded_states()
    for src, dest in WORKFLOWS.items():
        wf = inline_imports(load(src))
        native = python_to_workflow(wf, workflow_directory=SRC)
        pin_uuids(native, native["name"])
        report = []
        complete_tool_states(native, native["name"], recorded, report)
        for line in report:
            print(f"WARNING {dest}: {line}", file=sys.stderr)
        if report:
            print("  -> those steps keep the source state; Galaxy will fill their defaults and show the\n"
                  "     'tools have changed' notice. Refresh with tests/refresh_tool_states.py.", file=sys.stderr)
        with open(os.path.join(HERE, dest), "w") as f:
            json.dump(native, f, indent=4)
            f.write("\n")
        n_sub = sum(1 for s in native["steps"].values() if s["type"] == "subworkflow")
        print(f"{src} -> {dest}: {len(native['steps'])} steps ({n_sub} subworkflow)", file=sys.stderr)


if __name__ == "__main__":
    main()
