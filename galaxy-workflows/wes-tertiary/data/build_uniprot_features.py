#!/usr/bin/env python3
"""Protein features (UniProtKB/Swiss-Prot, CC BY 4.0) by Ensembl transcript, for the clinical report.

Input: the reviewed human UniProt entries as TSV from the UniProt REST API:
  https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=tsv&fields=accession,
  gene_primary,length,xref_ensembl,ft_domain,ft_region,ft_motif,ft_zn_fing,ft_dna_bind,ft_act_site,ft_binding,
  ft_repeat,ft_coiled
Feature positions are on UniProt's canonical isoform, so only the Ensembl transcripts that encode it (no isoform tag,
or the '-1' isoform) are kept; the report looks them up by its MANE transcript. Disordered regions predicted by
MobiDB-lite are left out.

Output columns: transcript (unversioned ENST)  gene  accession  start  end  type  note
Usage: build_uniprot_features.py uniprot.tsv uniprot_features.tsv.  Standard library only.
"""
import csv
import hashlib
import os
import re
import sys

FT = re.compile(r"(DOMAIN|REGION|MOTIF|ZN_FING|DNA_BIND|ACT_SITE|BINDING|REPEAT|COILED) (\d+)(?:\.\.(\d+))?;"
                r"(?: /note=\"([^\"]*)\";)?(?: /ligand=\"([^\"]*)\";)?")
TYPE = {"DOMAIN": "domain", "REGION": "region", "MOTIF": "motif", "ZN_FING": "zinc finger", "DNA_BIND": "DNA binding",
        "ACT_SITE": "active site", "BINDING": "binding site", "REPEAT": "repeat", "COILED": "coiled coil"}


def main(src, out):
    n = 0
    h = hashlib.sha256(open(src, "rb").read()).hexdigest()
    with open(src, encoding="utf-8", newline="") as f, open(out, "w", encoding="utf-8", newline="\n") as o:
        o.write(f"# source\t{os.path.basename(src)}\tsha256:{h}\n")
        o.write("transcript\tgene\taccession\tstart\tend\ttype\tnote\n")
        r = csv.reader(f, delimiter="\t")
        head = next(r)
        for t in r:
            row = dict(zip(head, t))
            acc = row["Entry"]
            tx = []
            for x in row.get("Ensembl", "").split(";"):
                m = re.match(r"\s*(ENST\d+)(?:\.\d+)?(?: \[([^\]]+)\])?", x)
                if m and (not m.group(2) or m.group(2) == f"{acc}-1"):
                    tx.append(m.group(1))
            if not tx:
                continue
            feats = []
            for col in head[4:]:
                for m in FT.finditer(row.get(col, "")):
                    note = m.group(4) or m.group(5) or ""
                    if note == "Disordered":
                        continue
                    s = int(m.group(2))
                    feats.append((s, int(m.group(3) or s), TYPE[m.group(1)], note))
            for e in dict.fromkeys(tx):
                for s, en, ty, note in feats:
                    o.write(f"{e}\t{row['Gene Names (primary)']}\t{acc}\t{s}\t{en}\t{ty}\t{note.replace(chr(9), ' ')}\n")
                    n += 1
    print(f"{n} feature rows", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
