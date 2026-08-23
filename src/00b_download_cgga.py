#!/usr/bin/env python3
"""Fetch CGGA mRNAseq_693 and mRNAseq_325 into data/raw/CGGA/.

Replication cohorts (CLAUDE.md 4):
  mRNAseq_693  -- MANDATORY.   GBM: 140 primary, 109 recurrent.
  mRNAseq_325  -- supportive.  GBM: 85 primary, 24 recurrent, 30 secondary.

For each cohort we want the expression matrix (read counts / RSEM gene table)
and the clinical table, which carries recurrence status and the covariates in
CLAUDE.md 4.3 (sequencing platform, IDH status, tumour purity).

Unlike GEO there is no derivable FTP path here, so links are scraped from the
CGGA download page. We do not hardcode or guess filenames. Run --dry-run first
and read the URLs before letting anything transfer.

TRANSPORT: if CGGA serves over plain HTTP, that is recorded in PROVENANCE.md.
The SHA256 still binds the file to exactly what we analysed and makes the
analysis reproducible -- but over HTTP it does not authenticate the source. Those
are different guarantees and the provenance record states which one we have.

Usage
-----
    python3 src/00b_download_cgga.py --dry-run
    python3 src/00b_download_cgga.py
"""

from __future__ import annotations

import argparse
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import (  # noqa: E402
    download_with_resume, extract_links, fetch_text, load_manifest, make_logger,
    remote_size, render_table, save_manifest, sha256_file, update_provenance, utc_now,
)

ACCESSION = "CGGA"
DOWNLOAD_PAGES = [
    "http://www.cgga.org.cn/download.jsp",
    "https://www.cgga.org.cn/download.jsp",
]

# Cohorts we need. A file qualifies if its name mentions the cohort AND looks
# like either an expression matrix or a clinical table.
COHORTS = ("mRNAseq_693", "mRNAseq_325")
EXPRESSION_HINTS = ("read_count", "readcount", "rsem", "fpkm", "count", "expression")
CLINICAL_HINTS = ("clinical",)

REPO = Path(__file__).resolve().parent.parent
DEST_DIR = REPO / "data" / "raw" / ACCESSION
PROVENANCE = REPO / "data" / "raw" / "PROVENANCE.md"
MANIFEST = DEST_DIR / ".download_manifest.json"
LOGFILE = REPO / "logs" / "00b_download_cgga.log"


def clean_name(url: str, raw: str) -> str:
    """Recover a real filename from a CGGA download URL.

    CGGA serves files through a servlet:
        download?file=download/20220620/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip&type=...
    so the last path segment carries the whole query string with it. The actual
    filename lives in the `file` query parameter.
    """
    query = urllib.parse.urlparse(url).query
    params = urllib.parse.parse_qs(query)
    if "file" in params and params["file"]:
        return params["file"][0].rsplit("/", 1)[-1]
    return raw.split("&")[0].split("?")[0]


def classify(name: str) -> tuple[str, str] | None:
    """Return (cohort, kind) if this filename is one we want, else None."""
    low = name.lower()
    for cohort in COHORTS:
        if cohort.lower().replace("_", "") in low.replace("_", "").replace(".", ""):
            if any(h in low for h in CLINICAL_HINTS):
                return cohort, "clinical"
            if any(h in low for h in EXPRESSION_HINTS):
                return cohort, "expression"
    return None


def discover(log) -> list[tuple[str, str, str, str]]:
    """Scrape the CGGA download page. Returns (name, url, cohort, kind)."""
    html, page = None, None
    for candidate in DOWNLOAD_PAGES:
        try:
            log(f"reading download page: {candidate}")
            html = fetch_text(candidate, timeout=60, retries=2, log=log)
            page = candidate
            break
        except Exception as e:
            log(f"  unreachable ({e})")
    if html is None:
        raise SystemExit("FATAL: could not reach the CGGA download page over HTTP or HTTPS.")

    log(f"page fetched over {page.split(':')[0].upper()}  ({len(html):,} bytes)")
    links = extract_links(html, page.rsplit("/", 1)[0])

    wanted: list[tuple[str, str, str, str]] = []
    seen: set[str] = set()
    for raw, url in links:
        name = clean_name(url, raw)
        hit = classify(name)
        if hit and name not in seen:
            seen.add(name)
            wanted.append((name, url, hit[0], hit[1]))

    log(f"page lists {len(links)} link(s); {len(wanted)} match the cohorts we need")
    return sorted(wanted, key=lambda x: (x[2], x[3], x[0]))


def report(files, log) -> None:
    if not files:
        log("NO MATCHING FILES FOUND — do not proceed. Inspect the page manually.")
        return
    log("")
    log(f"{'cohort':<14} {'kind':<11} {'transport':<9} {'bytes':>15}  file")
    log("-" * 100)
    for name, url, cohort, kind in files:
        size = remote_size(url, timeout=60)
        transport = "HTTPS" if url.startswith("https://") else "HTTP"
        shown = f"{size:,}" if size else "unreported"
        log(f"{cohort:<14} {kind:<11} {transport:<9} {shown:>15}  {name}")
        log(f"{'':>52}  {url}")
    log("")
    for cohort in COHORTS:
        kinds = {k for _, _, c, k in files if c == cohort}
        missing = {"expression", "clinical"} - kinds
        status = "complete" if not missing else f"MISSING {sorted(missing)}"
        log(f"  {cohort}: {status}")


def transport_note(files) -> str:
    """A PROVENANCE.md note describing what the checksums do and do not prove."""
    if any(u.startswith("http://") for _, u, _, _ in files):
        return (
            "> **Transport: plain HTTP.** CGGA served one or more of these files over\n"
            "> unencrypted HTTP. The SHA256 values below still bind each file to exactly\n"
            "> the bytes we analysed, so the analysis is reproducible and any later\n"
            "> alteration is detectable. They do **not** authenticate the source: an HTTP\n"
            "> transfer cannot prove the bytes came from CGGA unmodified. These are\n"
            "> different guarantees and only the first one is claimed here.\n"
        )
    return "> Transport: HTTPS, certificate verified.\n"


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
    log("00b_download_cgga.py starting")
    log(f"destination: {DEST_DIR}")

    files = discover(log)
    report(files, log)

    if args.dry_run:
        log("--dry-run: nothing downloaded")
        return 0
    if not files:
        return 1

    manifest = load_manifest(MANIFEST)
    ok, skipped, failed = 0, 0, 0

    for name, url, cohort, kind in files:
        dest = DEST_DIR / name
        log(f"{name}  [{cohort}/{kind}]")

        if dest.exists():
            recorded = manifest.get(name)
            digest = sha256_file(dest)
            if recorded and recorded["sha256"] == digest:
                log(f"  present, checksum matches manifest -- skip ({dest.stat().st_size:,} B)")
                skipped += 1
                continue
            if recorded and not args.force:
                log("  PRESENT BUT CHECKSUM MISMATCH — left untouched (data/raw is immutable).")
                log(f"    recorded: {recorded['sha256']}")
                log(f"    on disk:  {digest}")
                log("    re-run with --force to replace.")
                failed += 1
                continue
            if recorded:
                dest.unlink()
            else:
                log(f"  present, not in manifest -- hashing and recording")
                manifest[name] = {"file": name, "accession": f"CGGA {cohort} ({kind})",
                                  "url": url, "downloaded_utc": utc_now(),
                                  "bytes": dest.stat().st_size, "sha256": digest}
                save_manifest(MANIFEST, manifest)
                ok += 1
                continue

        try:
            log(f"  downloading {url}")
            n = download_with_resume(url, dest, timeout=args.timeout,
                                     retries=args.retries, log=log)
            digest = sha256_file(dest)
            manifest[name] = {"file": name, "accession": f"CGGA {cohort} ({kind})",
                              "url": url, "downloaded_utc": utc_now(),
                              "bytes": n, "sha256": digest}
            save_manifest(MANIFEST, manifest)
            log(f"  done: {n:,} B  sha256={digest}")
            ok += 1
        except Exception as e:
            log(f"  FAILED: {e}")
            failed += 1

    block = transport_note(files) + "\n" + render_table(list(manifest.values()))
    update_provenance(PROVENANCE, ACCESSION, block)
    log(f"provenance written: {PROVENANCE}")
    log(f"summary: {ok} downloaded/recorded, {skipped} skipped, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
