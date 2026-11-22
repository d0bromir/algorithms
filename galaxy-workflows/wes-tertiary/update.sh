#!/usr/bin/env bash
# Refresh the reference and knowledge files: build a new dated data release next to the current one, with only the
# sources that change re-downloaded, and (optionally) put it into Galaxy. Meant to run unattended (cron or the
# systemd timer in install/systemd/); everything stays on this server.
#
#   ./update.sh [options]
#
# 1. Takes the current release (ROOT/current, else the newest ROOT/release-*/ with a MANIFEST.tsv) and makes
#    ROOT/release-YYYY-MM-DD with hard links to all its files (no extra disk space).
# 2. Removes the files of the sources that change (FREQUENT below: ClinVar, GenCC, HPO, ClinGen, G2P/Orphanet,
#    UniProt, IPD-IMGT/HLA, IPD-KIR and the tables built from them) and runs "install.sh data" on the new release,
#    which downloads and rebuilds only what is missing. Version-pinned files (genome, Ensembl CDS, gnomAD, MANE,
#    AlphaMissense, REVEL, Collins, MaveDB, segmental duplications) stay linked; change them with install.sh.
# 3. Compares the refreshed files with the current release (content, not timestamps). Nothing changed: the new
#    release is removed. Otherwise it gets CHANGES.tsv (old and new date, size and source of each changed file),
#    becomes ROOT/current, and with --upload goes into a new Galaxy history ("WES reference data <release>"), where
#    unchanged files are copied from the previous release's history instead of uploaded again.
# 4. With --keep N, deletes all but the N newest releases (the data in Galaxy histories is not affected).
#
# Options:
#   --root DIR        directory of the releases (default /data/galaxy/wes-ref)
#   --regions BED     capture-kit regions, as for install.sh (needed to rebuild carrier_genes.tsv)
#   --upload          put a changed release into Galaxy (API key: $GALAXY_API_KEY or ~/.galaxy_admin_key)
#   --galaxy URL      Galaxy server (default http://localhost)
#   --keep N          keep only the N newest releases (default 0: keep all)
#   --force           keep the new release even when nothing changed
#   --dry-run         list what would be refreshed, change nothing
#
# Exit status: 0 (new release, or nothing changed, or another update still running), 1 (error; the partial release
# is removed and the current one is untouched). The log is ROOT/logs/update-YYYY-MM-DD.log and goes to stdout.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT=/data/galaxy/wes-ref
REGIONS=""
UPLOAD=0
GALAXY_URL="http://localhost"
KEEP=0
FORCE=0
DRY=0
# Sources that change between releases, and the files built from them (install.sh data rebuilds what is missing).
FREQUENT=(clinvar.GRCh38.vcf.gz clinvar_protein.tsv.gz gencc-submissions.tsv phenotype_to_genes.txt phenotype.hpoa
          ClinGen_gene_curation_list_GRCh38.tsv ClinGen_region_curation_list_GRCh38.tsv gene_disease.tsv
          vcep_rules.tsv erepo_variants.tsv uniprot_features.tsv.gz hla.dat kir.dat carrier_genes.tsv)

while (($#)); do
    case "$1" in
        --root) ROOT="$2"; shift 2 ;;
        --regions) REGIONS="$(realpath "$2")"; shift 2 ;;
        --upload) UPLOAD=1; shift ;;
        --galaxy) GALAXY_URL="$2"; shift 2 ;;
        --keep) KEEP="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0 ;;
        *) echo "update.sh: unknown argument '$1' (see --help)" >&2; exit 1 ;;
    esac
done

mkdir -p "$ROOT/logs"
LOG="$ROOT/logs/update-$(date +%F).log"
exec > >(tee -a "$LOG") 2>&1
echo "== update.sh $(date '+%F %T') on $(hostname)"

exec 9>"$ROOT/.update.lock"
flock -n 9 || { echo "another update is running: nothing done"; exit 0; }

current() {
    if [[ -L "$ROOT/current" && -f "$ROOT/current/MANIFEST.tsv" ]]; then readlink -f "$ROOT/current"; return; fi
    ls -1dt "$ROOT"/release-*/MANIFEST.tsv 2>/dev/null | head -1 | xargs -r dirname
}
PREV="$(current)"
[[ -n "$PREV" ]] || { echo "no release with a MANIFEST.tsv under $ROOT: run install.sh data first"; exit 1; }
NEW="$ROOT/release-$(date +%F)"
n=2; while [[ -e "$NEW" ]]; do NEW="$ROOT/release-$(date +%F)-$n"; n=$((n + 1)); done
echo "current release: $PREV"
echo "new release:     $NEW"

if [[ -z "$REGIONS" && -s "$PREV/regions.bed" ]]; then REGIONS="$PREV/regions.bed"; fi
if ((DRY)); then
    for f in "${FREQUENT[@]}"; do
        [[ -e "$PREV/$f" ]] || continue
        if [[ $f == carrier_genes.tsv && -z "$REGIONS" ]]; then echo "  would keep $f (no --regions)"
        else echo "  would refresh $f"; fi
    done
    exit 0
fi

fail() { echo "update.sh: $*; removing $NEW"; rm -rf "$NEW"; exit 1; }
trap 'fail "interrupted or failed at line $LINENO"' ERR

# install.sh writes new files (temporary name, then rename) and leaves present ones alone, so the hard-linked files
# of the current release are never modified through the new one.
cp -al "$PREV" "$NEW"
rm -f "$NEW/CHANGES.tsv"
for f in "${FREQUENT[@]}"; do
    if [[ $f == carrier_genes.tsv && -z "$REGIONS" ]]; then
        echo "  no --regions and no regions.bed: carrier_genes.tsv kept from $(basename "$PREV")"; continue
    fi
    rm -f "$NEW/$f"
done
args=(--data-dir "$NEW")
[[ -n "$REGIONS" ]] && args+=(--regions "$REGIONS")
"$HERE/install.sh" "${args[@]}" data || fail "install.sh data failed"

content() {  # hash of the data lines: comment lines (#) carry download dates and source hashes, and gzip headers a
             # timestamp, so neither counts as a change
    case "$1" in
        *.gz) gzip -dc "$1" | grep -v '^#' | sha256sum | cut -d' ' -f1 ;;
        *) grep -v '^#' "$1" | sha256sum | cut -d' ' -f1 ;;
    esac
}
row() { awk -F'\t' -v f="$2" '$1 == f {print $4"\t"$5"\t"$3}' "$1/MANIFEST.tsv"; }
changed=()
for f in "${FREQUENT[@]}"; do
    [[ -e "$NEW/$f" ]] || continue
    if [[ ! -e "$PREV/$f" ]] || [[ "$(content "$PREV/$f")" != "$(content "$NEW/$f")" ]]; then changed+=("$f"); fi
done
trap - ERR
# install.sh keeps a file's date when its size and source are unchanged; a changed file is dated today
for f in "${changed[@]}"; do
    awk -F'\t' -v OFS='\t' -v f="$f" -v d="$(date +%F)" '$1 == f {$4 = d} {print}' "$NEW/MANIFEST.tsv" > "$NEW/MANIFEST.tsv.tmp"
    mv "$NEW/MANIFEST.tsv.tmp" "$NEW/MANIFEST.tsv"
done

if ((${#changed[@]} == 0 && ! FORCE)); then
    echo "nothing changed since $(basename "$PREV"): new release removed"
    rm -rf "$NEW"
    exit 0
fi
{
    printf 'file\tprevious_date\tprevious_bytes\tprevious_source\tnew_date\tnew_bytes\tnew_source\n'
    for f in "${changed[@]}"; do
        p="$(row "$PREV" "$f")"; [[ -n "$p" ]] || p=$'-\t-\t-'
        printf '%s\t%s\t%s\n' "$f" "$p" "$(row "$NEW" "$f")"
    done
} > "$NEW/CHANGES.tsv"
echo "changed: ${changed[*]:-none (--force)}"
column -t -s $'\t' "$NEW/CHANGES.tsv" | sed 's/^/  /'
ln -sfn "$(basename "$NEW")" "$ROOT/current"
echo "current -> $(basename "$NEW")"

if ((UPLOAD)); then
    "$HERE/install.sh" --galaxy "$GALAXY_URL" --data-dir "$NEW" \
        --previous-history "WES reference data $(basename "$PREV")" upload ||
        { echo "update.sh: upload failed (the release is kept; re-run: install.sh --data-dir $NEW upload)"; exit 1; }
fi

if ((KEEP > 0)); then
    mapfile -t old < <(ls -1dt "$ROOT"/release-*/ | tail -n +$((KEEP + 1)))
    for d in "${old[@]}"; do
        [[ "$(readlink -f "$d")" == "$(readlink -f "$NEW")" ]] && continue
        echo "removing old release $d"; rm -rf "$d"
    done
fi
echo "== done $(date '+%F %T')"
