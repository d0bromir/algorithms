#!/usr/bin/env python3
"""Galaxy side of install.sh (BioBlend and Ephemeris). The admin API key comes from $GALAXY_API_KEY.

    galaxy_setup.py check TOOLS_YML TOOL_ID...      what is installed: Tool Shed repositories, local tools, workflows
    galaxy_setup.py tools TOOLS_YML                 install the Tool Shed repositories (Ephemeris shed-tools)
    galaxy_setup.py local-deps TOOL_ID...           install the Conda dependencies of local tools that lack them
    galaxy_setup.py workflows [--update] GA...      import workflows (or update those with the same name)
    galaxy_setup.py upload DATA_DIR HISTORY SNPEFF  upload the MANIFEST.tsv files to a new history, fetch SnpEff db

Galaxy URL: $GALAXY_URL (default http://localhost).
"""
import csv
import json
import os
import sys
import time

import yaml
from bioblend.galaxy import GalaxyInstance

URL = os.environ.get("GALAXY_URL", "http://localhost")
KEY = os.environ.get("GALAXY_API_KEY") or sys.exit("GALAXY_API_KEY is not set")
gi = GalaxyInstance(URL, key=KEY)


def installed_repos():
    return {(r["name"], r["owner"]) for r in gi.toolshed.get_repositories()
            if r["status"] == "Installed" and not r["deleted"]}


def check(tools_yml, *tool_ids):
    want = [(t["name"], t["owner"]) for t in yaml.safe_load(open(tools_yml))["tools"]]
    have = installed_repos()
    missing = [f"{o}/{n}" for n, o in want if (n, o) not in have]
    print(f"  Tool Shed repositories: {len(want) - len(missing)} of {len(want)} installed"
          + (f"; missing: {', '.join(missing)}" if missing else ""))
    for tid in tool_ids:
        try:
            t = gi.tools.show_tool(tid)
        except Exception:  # noqa: BLE001 - not loaded
            print(f"  local tool {tid}: not loaded by Galaxy")
            continue
        reqs = gi.tools.requirements(tid)
        bad = [f"{r['name']} {r['version']}" for r in reqs if r.get("dependency_type") in (None, "null")]
        print(f"  local tool {tid} {t['version']}: dependencies "
              + (f"MISSING ({', '.join(bad)})" if bad else f"resolved ({len(reqs)})"))
    names = {w["name"] for w in gi.workflows.get_workflows()}
    for n in ("WES tertiary analysis from VCF (v2)", "WES tertiary analysis from FASTQ (v1)"):
        print(f"  workflow {n}: {'imported' if n in names else 'not imported'}")


def tools(tools_yml):
    """Ephemeris in-process: the key stays out of the process list."""
    from ephemeris import shed_tools
    want = yaml.safe_load(open(tools_yml))["tools"]
    have = installed_repos()
    todo = [t for t in want if (t["name"], t["owner"]) not in have]
    if not todo:
        print(f"  all {len(want)} repositories installed")
        return
    print(f"  installing {len(todo)}: {', '.join(t['owner'] + '/' + t['name'] for t in todo)}")
    sys.argv = ["shed-tools", "install", "-g", URL, "-a", KEY, "-t", tools_yml]
    try:
        shed_tools.main()
    except SystemExit as e:
        if e.code:
            raise
    left = [t["name"] for t in todo if (t["name"], t["owner"]) not in installed_repos()]
    if left:
        sys.exit(f"  not installed: {', '.join(left)}")


def local_deps(*tool_ids):
    for tid in tool_ids:
        try:
            gi.tools.show_tool(tid)
        except Exception:  # noqa: BLE001 - not loaded
            print(f"  {tid}: not loaded by Galaxy yet: restart Galaxy, then re-run the 'local-tools' stage")
            sys.exit(1)
        reqs = gi.tools.requirements(tid)
        if reqs and all(r.get("dependency_type") not in (None, "null") for r in reqs):
            print(f"  {tid}: resolved ({reqs[0].get('environment_path', '')})")
            continue
        print(f"  {tid}: installing {', '.join(r['name'] + ' ' + r['version'] for r in reqs)} ...", flush=True)
        t0 = time.time()
        gi.tools.install_dependencies(tid)
        reqs = gi.tools.requirements(tid)
        bad = [r["name"] for r in reqs if r.get("dependency_type") in (None, "null")]
        print(f"  {tid}: {'FAILED for ' + ', '.join(bad) if bad else 'installed'} in {time.time() - t0:.0f} s")
        if bad:
            sys.exit(1)
    print("  To reclaim the Conda package cache (it grows to 10+ GB), run as the Galaxy user: conda clean -a -y")


def count(m):
    if isinstance(m, dict):
        return sum(count(v) for v in m.values())
    if isinstance(m, list):
        return sum(count(v) for v in m)
    return 1 if m else 0


def workflows(*args):
    update = "--update" in args
    existing = {w["name"]: w for w in gi.workflows.get_workflows()}
    for path in (a for a in args if a != "--update"):
        native = json.load(open(path))
        name = native["name"]
        if name in existing and not update:
            print(f"  {name}: present (use --update for a new version from {os.path.basename(path)})")
            continue
        if name in existing:
            wid = existing[name]["id"]
            gi.workflows._put(id=wid, payload={"workflow": native})  # new version of the stored workflow
            what = "updated"
        else:
            wid = gi.workflows.import_workflow_dict(native)["id"]
            what = "imported"
        ed = gi.workflows._get(id=f"{wid}/download", params={"style": "editor"})
        print(f"  {name}: {what} (release {native.get('release')}), "
              f"upgrade messages {count(ed.get('upgrade_messages'))}")


def wait(dataset_ids, what):
    t0 = time.time()
    for did in dataset_ids:
        while True:
            st = gi.datasets.show_dataset(did)["state"]
            if st in ("ok", "error", "failed_metadata", "discarded", "paused"):
                break
            time.sleep(10)
        if st != "ok":
            sys.exit(f"  {what}: dataset {did} ended in state {st}")
    print(f"  {what}: ok ({time.time() - t0:.0f} s)")


def upload(data_dir, history, snpeff_db):
    rows = list(csv.DictReader(open(os.path.join(data_dir, "MANIFEST.tsv")), delimiter="\t"))
    h = gi.histories.create_history(history)
    print(f"  history '{history}' ({h['id']})")
    ids = []
    for r in rows:
        p = os.path.join(data_dir, r["file"])
        out = gi.tools.upload_file(p, h["id"], file_type=r["galaxy_datatype"], file_name=r["file"])  # tus upload
        ids += [o["id"] for o in out["outputs"]]
        print(f"  uploading {r['file']} ({int(r['bytes']) / 1024**2:,.0f} MB) as {r['galaxy_datatype']}", flush=True)
    dl = [t["id"] for t in gi.tools.get_tools() if "/snpeff/snpEff_download/" in t["id"]]
    if dl:
        tid = sorted(dl, key=lambda x: [int(p) if p.isdigit() else p for p in x.split("/")[-1].replace("+", ".").split(".")])[-1]
        out = gi.tools.run_tool(h["id"], tid, {"genome_version": snpeff_db})
        ids += [o["id"] for o in out["outputs"]]
        print(f"  SnpEff database {snpeff_db} (tool {tid.split('/repos/')[-1]})")
    else:
        print("  SnpEff download tool not installed: run the 'tools' stage, then the SnpEff download tool")
    wait(ids, f"{len(ids)} datasets")


if __name__ == "__main__":
    cmds = {"check": check, "tools": tools, "local-deps": local_deps, "workflows": workflows, "upload": upload}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        sys.exit(__doc__)
    cmds[sys.argv[1]](*sys.argv[2:])
