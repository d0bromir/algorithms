#!/usr/bin/env python3
"""Record the complete tool states Galaxy computes for the WES workflows in src/galaxy-tool-states.json.

The gxformat2 sources only set the parameters that matter; Galaxy fills every other parameter with its
default when it loads a step and then flags the workflow ("may contain tools which have changed").
build.py merges the states recorded here into the .ga files, so imports open without that notice.

Refresh after changing a source or a tool version:
    1. python build.py                      (build the .ga files)
    2. import both .ga files into a Galaxy that has the tools of tools.yml installed
    3. python tests/refresh_tool_states.py --galaxy http://192.168.191.12   (key: $GALAXY_API_KEY
       or ~/.galaxy_admin_key)
    4. python build.py                      (now with the full states)
The workflows are looked up by name among the API user's workflows (most recently updated first).
Only Python's standard library is needed.
"""
import argparse
import json
import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
WF_DIR = os.path.dirname(HERE)
STATES = os.path.join(WF_DIR, "src", "galaxy-tool-states.json")
GA_FILES = ["wes-tertiary-from-vcf.ga", "wes-tertiary-from-fastq.ga"]


def decode(tool_state):
    """Galaxy's editor state (top-level values JSON-encoded) -> plain nested dict."""
    if isinstance(tool_state, str):
        tool_state = json.loads(tool_state)
    return {k: (json.loads(v) if isinstance(v, str) and not k.startswith("__") else v)
            for k, v in tool_state.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxy", required=True)
    ap.add_argument("--api-key-file", default=os.path.expanduser("~/.galaxy_admin_key"))
    args = ap.parse_args()
    key = os.environ.get("GALAXY_API_KEY") or "".join(open(args.api_key_file).read().split())

    def get(path):
        req = urllib.request.Request(args.galaxy.rstrip("/") + path, headers={"x-api-key": key})
        with urllib.request.urlopen(req) as r:
            return json.loads(r.read())

    by_name = {}
    for w in get("/api/workflows?show_deleted=false"):
        by_name.setdefault(w["name"], w)  # newest first

    editor = {}
    states = {}
    for fname in GA_FILES:
        native = json.load(open(os.path.join(WF_DIR, fname)))
        name = native["name"]
        if name not in by_name:
            raise SystemExit(f"'{name}' is not imported in {args.galaxy}")
        ed = get(f"/api/workflows/{by_name[name]['id']}/download?style=editor")
        editor[name] = {s["label"]: s for s in ed["steps"].values() if s.get("type") == "tool"}

    def record(native, scope, source):
        for s in native["steps"].values():
            if s.get("type") == "tool":
                es = editor[source][s["label"]]
                states[f"{scope}/{s['label']}"] = {
                    "tool_id": es["content_id"],
                    "tool_version": s["tool_version"],
                    "tool_state": decode(es["tool_state"]),
                }
            elif s.get("type") == "subworkflow":
                # The embedded workflow is also imported on its own; Galaxy's states for it are the same.
                record(s["subworkflow"], f"{scope}/{s['label']}", s["subworkflow"]["name"])

    for fname in GA_FILES:
        native = json.load(open(os.path.join(WF_DIR, fname)))
        record(native, native["name"], native["name"])

    with open(STATES, "w") as f:
        json.dump(states, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"wrote {len(states)} step states to {os.path.relpath(STATES, WF_DIR)}")


if __name__ == "__main__":
    main()
