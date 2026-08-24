#!/usr/bin/env python3
"""Fetch hg38 gene coordinates and cytobands -> data/raw/ucsc_hg38/.

infercnvpy needs `chromosome`, `start`, `end` on adata.var. We additionally need
the p/q ARM of each gene, because disjoint_set_S is arm-aware: chr9**p** only, with
9q genes remaining eligible for Stage B (CLAUDE.md §3.4).

Two open, direct-download files from UCSC (no account, no DUA):
  refFlat.txt.gz   gene symbol -> chrom, txStart, txEnd
  cytoBand.txt.gz  cytogenetic bands, including the `acen` centromere bands that
                   define the p/q boundary

Downloads and records only. No parsing into the pipeline here.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import (  # noqa: E402
    NETWORK_ERRORS, download_with_resume, load_manifest, make_logger, render_table,
    save_manifest, sha256_file, update_provenance, utc_now,
)

ACCESSION = "ucsc_hg38"
BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
FILES = {"refFlat.txt.gz": f"{BASE}/refFlat.txt.gz",
         "cytoBand.txt.gz": f"{BASE}/cytoBand.txt.gz"}

REPO = Path(__file__).resolve().parent.parent
DEST = REPO / "data" / "raw" / ACCESSION
PROVENANCE = REPO / "data" / "raw" / "PROVENANCE.md"
MANIFEST = DEST / ".download_manifest.json"
LOGFILE = REPO / "logs" / "00g_download_gene_annotation.log"

NOTE = ("> **UCSC hg38 annotation.** Open access, no account or DUA. `refFlat`\n"
        "> supplies gene coordinates for inferCNV windows; `cytoBand` supplies the\n"
        "> `acen` centromere bands that define the p/q arm boundary, which\n"
        "> `disjoint_set_S` requires because chr9**p** is an arm, not a chromosome.\n")


def main() -> int:
    DEST.mkdir(parents=True, exist_ok=True)
    log = make_logger(LOGFILE)
    log("=" * 72)
    log("00g_download_gene_annotation.py starting")
    manifest = load_manifest(MANIFEST)
    failed = 0

    for name, url in FILES.items():
        dest = DEST / name
        if dest.exists():
            digest = sha256_file(dest)
            rec = manifest.get(name)
            if rec and rec["sha256"] == digest:
                log(f"{name}: present, checksum matches -- skip")
                continue
        try:
            log(f"{name}: downloading {url}")
            n = download_with_resume(url, dest, log=log)
            digest = sha256_file(dest)
            manifest[name] = {"file": name, "accession": "UCSC hg38",
                              "url": url, "downloaded_utc": utc_now(),
                              "bytes": n, "sha256": digest}
            save_manifest(MANIFEST, manifest)
            log(f"  done: {n:,} B  sha256={digest}")
        except (NETWORK_ERRORS + (IOError,)) as e:
            log(f"  FAILED: {e}")
            failed += 1

    update_provenance(PROVENANCE, ACCESSION,
                      NOTE + "\n" + render_table(list(manifest.values())))
    log(f"provenance written; {failed} failure(s)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
