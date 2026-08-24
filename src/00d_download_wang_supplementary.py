#!/usr/bin/env python3
"""Fetch Wang et al. 2022 supplementary files into data/raw/wang2022_supplementary/.

    Wang L, Jung J, Babikir H, et al.
    "A single-cell atlas of glioblastoma evolution under therapy reveals
     cell-intrinsic and cell-extrinsic therapeutic targets."
    Nature Cancer 3(12):1534-1552, 2022.
    PMID 36539501 | DOI 10.1038/s43018-022-00475-x | PMC9767870 | licence CC BY

Why: per-patient IDH status is NOT in the GEO deposit. GSE174554's only sample
characteristics are progression, diagnosis, age, gender, pair# and tissue
(verified against GSE174554_family.soft.gz). Any IDH claim must come from the
paper's supplementary tables — see CLAUDE.md §4.

Route: www.nature.com serves a bot stub (~3 KB, no ESM links) to non-browser
clients, so we use the Europe PMC REST supplementaryFiles endpoint, which returns
the complete supplementary bundle as one zip. The zip is stored byte-exact in
data/raw/ with a checksum; extraction goes to data/interim/ (CLAUDE.md §7.10).

This script downloads and records. It does not parse any table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import (  # noqa: E402
    NETWORK_ERRORS, download_with_resume, load_manifest, make_logger, remote_size,
    render_table, save_manifest, sha256_file, update_provenance, utc_now,
)

ACCESSION = "wang2022_supplementary"
PMCID = "PMC9767870"
URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
FILENAME = f"{PMCID}_supplementaryFiles.zip"

REPO = Path(__file__).resolve().parent.parent
DEST_DIR = REPO / "data" / "raw" / ACCESSION
PROVENANCE = REPO / "data" / "raw" / "PROVENANCE.md"
MANIFEST = DEST_DIR / ".download_manifest.json"
LOGFILE = REPO / "logs" / "00d_download_wang_supplementary.log"

CITATION_NOTE = (
    "> **Wang et al. 2022 supplementary bundle.** PMID 36539501, DOI\n"
    "> 10.1038/s43018-022-00475-x, PMC9767870, **licence CC BY**. Retrieved from the\n"
    "> Europe PMC REST `supplementaryFiles` endpoint over HTTPS because\n"
    "> `www.nature.com` returns a ~3 KB bot stub with no supplementary links to\n"
    "> non-browser clients. The archive is stored exactly as served; extraction goes\n"
    "> to `data/interim/` (CLAUDE.md §7.10). **Source of per-patient IDH status,\n"
    "> which the GEO deposit does not carry.**\n"
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--retries", type=int, default=5)
    args = ap.parse_args()

    DEST_DIR.mkdir(parents=True, exist_ok=True)
    log = make_logger(LOGFILE)
    log("=" * 72)
    log("00d_download_wang_supplementary.py starting")
    log(f"  {URL}")

    size = remote_size(URL, timeout=args.timeout)
    log(f"  remote size: {size:,} B" if size else "  remote size: unreported")
    if args.dry_run:
        log("--dry-run: nothing downloaded")
        return 0

    manifest = load_manifest(MANIFEST)
    dest = DEST_DIR / FILENAME

    if dest.exists():
        digest = sha256_file(dest)
        rec = manifest.get(FILENAME)
        if rec and rec["sha256"] == digest and not args.force:
            log(f"present, checksum matches manifest -- skip ({dest.stat().st_size:,} B)")
            update_provenance(PROVENANCE, ACCESSION,
                              CITATION_NOTE + "\n" + render_table(list(manifest.values())))
            return 0
        if rec and not args.force:
            log("PRESENT BUT CHECKSUM MISMATCH — left untouched (data/raw is immutable).")
            log("  re-run with --force to replace.")
            return 1
        dest.unlink()

    try:
        n = download_with_resume(URL, dest, timeout=args.timeout,
                                 retries=args.retries, log=log)
    except (NETWORK_ERRORS + (IOError,)) as e:
        log(f"FAILED: {e}")
        return 1

    digest = sha256_file(dest)
    manifest[FILENAME] = {
        "file": FILENAME,
        "accession": f"Wang et al. 2022 ({PMCID}, CC BY)",
        "url": URL, "downloaded_utc": utc_now(), "bytes": n, "sha256": digest,
    }
    save_manifest(MANIFEST, manifest)
    log(f"done: {n:,} B  sha256={digest}")

    update_provenance(PROVENANCE, ACCESSION,
                      CITATION_NOTE + "\n" + render_table(list(manifest.values())))
    log(f"provenance written: {PROVENANCE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
