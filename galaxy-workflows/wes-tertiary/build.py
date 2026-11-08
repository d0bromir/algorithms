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


def main():
    for src, dest in WORKFLOWS.items():
        wf = inline_imports(load(src))
        native = python_to_workflow(wf, workflow_directory=SRC)
        pin_uuids(native, native["name"])
        with open(os.path.join(HERE, dest), "w") as f:
            json.dump(native, f, indent=4)
            f.write("\n")
        n_sub = sum(1 for s in native["steps"].values() if s["type"] == "subworkflow")
        print(f"{src} -> {dest}: {len(native['steps'])} steps ({n_sub} subworkflow)", file=sys.stderr)


if __name__ == "__main__":
    main()
