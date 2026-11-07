#!/usr/bin/env python3
"""Build the importable Galaxy workflows (.ga) from the gxformat2 sources in ./src.

    pip install gxformat2
    python build.py

* ``run: {"@import": "<file>.gxwf.yml"}`` in a step is replaced by the content of that file, so the
  FASTQ workflow carries the VCF workflow as an embedded (self-contained) subworkflow.
* UUIDs are derived from the workflow/step labels (uuid5), so rebuilding an unchanged source produces
  an identical .ga file.
"""
import copy
import json
import os
import sys
import uuid

from gxformat2.converter import python_to_workflow
from gxformat2.yaml import ordered_load

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


def pin_uuids(native, scope):
    native["uuid"] = str(uuid.uuid5(NAMESPACE, scope))
    for key, step in native["steps"].items():
        step_scope = f"{scope}/{step.get('label') or key}"
        if "uuid" in step:
            step["uuid"] = str(uuid.uuid5(NAMESPACE, step_scope))
        for i, wo in enumerate(step.get("workflow_outputs") or []):
            if "uuid" in wo:
                wo["uuid"] = str(uuid.uuid5(NAMESPACE, f"{step_scope}/output/{wo.get('label') or i}"))
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
