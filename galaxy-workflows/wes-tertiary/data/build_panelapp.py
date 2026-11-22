#!/usr/bin/env python3
"""Snapshot of the gene panels of PanelApp (Genomics England) and PanelApp Australia, for offline use.

    build_panelapp.py panelapp_panels.tsv [--source ge,au] [--signed-off]

Reads each public panel once over HTTPS (only panel numbers are requested; no sample data leaves the server) and
writes one row per panel gene:
    source  panel_id  panel  version  version_created  signed_off  disorders  gene  hgnc_id  rating  moi  types
rating: green (confidence 3, diagnostic grade), amber (2), red (1); types: PanelApp's panel types (Super Panel,
Component Of Super Panel, GMS signed-off ...), which mark the broad superpanels. Kept are the panel, its version and
test codes (relevant disorders, e.g. R169), the gene and its rating and mode of inheritance; the phenotype text is left out
(part of it derives from OMIM, whose licence PanelApp's terms of use pass on). PanelApp's terms: non-commercial use,
with acknowledgement. Short tandem repeats and regions are not kept (the report works on genes).

Rows are sorted, so a rebuild from unchanged panels gives the same file (update.sh compares contents).
Standard library only.
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.request

API = {"ge": ("Genomics England PanelApp", "https://panelapp.genomicsengland.co.uk/api/v1"),
       "au": ("PanelApp Australia", "https://panelapp-aus.org/api/v1")}
RATING = {"3": "green", "2": "amber", "1": "red", "4": "green"}


DELAY = 1.0  # seconds between requests: PanelApp answers 429 (Too Many Requests) to faster clients


def get(url, tries=8):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "wes-tertiary"})
            with urllib.request.urlopen(req, timeout=60) as r:
                time.sleep(DELAY)
                return json.load(r)
        except urllib.error.HTTPError as e:  # rate limit: wait as told (Retry-After), else back off
            if i == tries - 1:
                raise SystemExit(f"{url}: {e}")
            time.sleep(float(e.headers.get("Retry-After") or 0) or min(2 ** (i + 2), 120))
        except Exception as e:  # transient network errors
            if i == tries - 1:
                raise SystemExit(f"{url}: {e}")
            time.sleep(min(2 ** (i + 2), 120))


def pages(url):
    while url:
        d = get(url)
        yield from d.get("results", [])
        url = d.get("next")


def clean(x):
    return " ".join(str(x or "").split()).replace("\t", " ")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out")
    ap.add_argument("--source", default="ge,au", help="ge, au or both (default)")
    ap.add_argument("--signed-off", action="store_true", help="Genomics England: only the signed-off panels")
    a = ap.parse_args(argv)
    rows, n_panels = [], 0
    for src in a.source.split(","):
        name, base = API[src]
        signed = {p["id"]: p.get("version") for p in pages(f"{base}/panels/signedoff/")} if src == "ge" else {}
        panels = list(pages(f"{base}/panels/"))
        for p in panels:
            if a.signed_off and src == "ge" and p["id"] not in signed:
                continue
            d = get(f"{base}/panels/{p['id']}/")
            n_panels += 1
            for g in d.get("genes", []):
                gd = g.get("gene_data") or {}
                rows.append([src, str(d["id"]), clean(d["name"]), clean(d.get("version")),
                             clean(d.get("version_created", ""))[:10],
                             clean(signed.get(d["id"], "")) if src == "ge" else "",
                             ";".join(clean(x) for x in d.get("relevant_disorders") or []),
                             clean(gd.get("gene_symbol") or g.get("entity_name")), clean(gd.get("hgnc_id")),
                             RATING.get(str(g.get("confidence_level")), "red"),
                             clean(g.get("mode_of_inheritance")),
                             ";".join(sorted(clean(t.get("name")) for t in d.get("types") or [] if t.get("name")))])
        print(f"{name}: {len(panels)} panels", file=sys.stderr)
    rows.sort()
    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        f.write("# source\t" + ", ".join(API[s][0] for s in a.source.split(",")) +
                f"\t{time.strftime('%Y-%m-%d')}\tnon-commercial use, with acknowledgement (PanelApp terms of use)\n")
        f.write("source\tpanel_id\tpanel\tversion\tversion_created\tsigned_off\tdisorders\tgene\thgnc_id\trating\tmoi"
                "\ttypes\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{n_panels} panels, {len(rows)} panel genes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
