#!/usr/bin/env python3
"""Fetch GSE174554 supplementary files from GEO into data/raw/GSE174554/.

Discovery cohort: Wang et al. 2022, Nature Cancer. Matched primary/recurrent
GBM, single-nucleus RNA-seq. Open access, no account, no DUA (CLAUDE.md 4).

This script downloads bytes and records provenance. It does not parse,
decompress, or load anything. data/raw/ is immutable (CLAUDE.md 7.4): a file
already present with a matching recorded checksum is skipped, and a file
present with a MISMATCHING checksum is reported and left alone unless --force.

Usage
-----
    python3 src/00_download.py --dry-run     # list what GEO offers, fetch nothing
    python3 src/00_download.py               # download + record provenance
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _download_utils import (  # noqa: E402
    NETWORK_ERRORS, download_with_resume, extract_links, fetch_text, load_manifest,
    make_logger, remote_size, render_table, save_manifest, sha256_file,
    update_provenance, utc_now,
)

ACCESSION = "GSE174554"

# The authors' own malignant vs non-malignant annotation. CLAUDE.md 4 makes this
# mandatory -- "Use this. Do not roll our own classifier." If GEO does not list
# it, that is a protocol-level problem and the run stops rather than continuing
# with a partial download that looks complete.
# NOTE: CLAUDE.md 4 names this "GSE174554_Tumor_normal_metadata.txt"; GEO actually
# serves it gzipped. We match on the stem so either form satisfies the check.
REQUIRED_FILES = {"GSE174554_Tumor_normal_metadata.txt"}

# Server-side index pages, not deposit content.
NOT_DATA = {"index.html", "filelist.txt"}


def _stem(name: str) -> str:
    """Filename with a trailing .gz removed, for required-file matching."""
    return name[:-3] if name.endswith(".gz") else name

REPO = Path(__file__).resolve().parent.parent
DEST_DIR = REPO / "data" / "raw" / ACCESSION
PROVENANCE = REPO / "data" / "raw" / "PROVENANCE.md"
MANIFEST = DEST_DIR / ".download_manifest.json"
LOGFILE = REPO / "logs" / "00_download.log"


# GEO publishes a series under several sibling directories. suppl/ holds the
# authors' uploaded matrices; soft/ and matrix/ hold GEO's own per-sample records.
# The SOFT family file is the authoritative source for per-GSM characteristics --
# organism, platform, library strategy, and the sample-level fields we need to
# resolve patient and timepoint (CLAUDE.md 10.1).
SUBDIRS = ("suppl", "soft", "matrix")


def series_index_url(accession: str, subdir: str) -> str:
    """GEO nests series by accession prefix: GSE174554 -> GSE174nnn."""
    stem = accession[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/{subdir}/"


def discover(log) -> list[tuple[str, str]]:
    """Read the file list from every GEO subdirectory. Filenames are never guessed."""
    files: list[tuple[str, str]] = []
    failed_subdirs: list[str] = []
    for subdir in SUBDIRS:
        index = series_index_url(ACCESSION, subdir)
        log(f"reading {subdir}/ index: {index}")
        try:
            found = extract_links(fetch_text(index, log=log), index)
        except NETWORK_ERRORS as e:
            log(f"  {subdir}/ unavailable ({e}) -- skipping")
            failed_subdirs.append(subdir)
            continue
        # Drop server index pages and absolute links off-site (GEO footers).
        found = [(n, u) for n, u in found
                 if n not in NOT_DATA and u.startswith(index)]
        log(f"  {subdir}/ lists {len(found)} file(s)")
        files.extend(found)

    # A subdirectory silently going missing would quietly shrink the download and
    # still report success. suppl/ carries the required annotation, so its absence
    # is fatal; the others are reported loudly.
    if "suppl" in failed_subdirs:
        raise SystemExit("FATAL: GEO suppl/ index unreachable -- refusing to "
                         "continue with a partial file list.")
    if failed_subdirs:
        log(f"WARNING: unreachable subdirectories skipped: {failed_subdirs}")
    log(f"GEO lists {len(files)} file(s) across "
        f"{len(SUBDIRS) - len(failed_subdirs)}/{len(SUBDIRS)} subdirectories")
    for name, url in sorted(files):
        size = remote_size(url)
        log(f"  - {name}" + (f"  ({size:,} B)" if size else "  (size unreported)"))
    missing = REQUIRED_FILES - {_stem(n) for n, _ in files}
    if missing:
        raise SystemExit(
            f"FATAL: required file(s) not listed by GEO: {sorted(missing)}\n"
            f"CLAUDE.md 4 requires the authors' tumour/normal annotation. Stop and "
            f"check the accession before proceeding."
        )
    return sorted(files)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="list the files GEO offers and exit without downloading")
    ap.add_argument("--force", action="store_true",
                    help="re-download files whose checksum does not match the manifest")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--retries", type=int, default=5)
    args = ap.parse_args()

    DEST_DIR.mkdir(parents=True, exist_ok=True)
    log = make_logger(LOGFILE)
    log("=" * 72)
    log(f"00_download.py starting -- accession {ACCESSION}")
    log(f"destination: {DEST_DIR}")

    files = discover(log)
    if args.dry_run:
        log("--dry-run: nothing downloaded")
        return 0

    manifest = load_manifest(MANIFEST)
    ok, skipped, failed = 0, 0, 0

    for name, url in files:
        dest = DEST_DIR / name
        log(f"{name}")

        if dest.exists():
            recorded = manifest.get(name)
            digest = sha256_file(dest)
            if recorded and recorded["sha256"] == digest:
                log(f"  present, checksum matches manifest -- skip ({dest.stat().st_size:,} B)")
                skipped += 1
                continue
            if recorded:
                log(f"  PRESENT BUT CHECKSUM MISMATCH")
                log(f"    recorded: {recorded['sha256']}")
                log(f"    on disk:  {digest}")
                if not args.force:
                    log("    left untouched (data/raw is immutable). Re-run with --force to replace.")
                    failed += 1
                    continue
                log("    --force given: re-downloading")
                dest.unlink()
            else:
                # Present but never recorded, e.g. an interrupted earlier run
                # that renamed before writing the manifest. Adopt and record it.
                log(f"  present, not in manifest -- hashing and recording ({dest.stat().st_size:,} B)")
                manifest[name] = {"file": name, "accession": ACCESSION, "url": url,
                                  "downloaded_utc": utc_now(), "bytes": dest.stat().st_size,
                                  "sha256": digest}
                save_manifest(MANIFEST, manifest)
                ok += 1
                continue

        try:
            log(f"  downloading {url}")
            n = download_with_resume(url, dest, timeout=args.timeout,
                                     retries=args.retries, log=log)
            digest = sha256_file(dest)
            manifest[name] = {"file": name, "accession": ACCESSION, "url": url,
                              "downloaded_utc": utc_now(), "bytes": n, "sha256": digest}
            save_manifest(MANIFEST, manifest)
            log(f"  done: {n:,} B  sha256={digest}")
            ok += 1
        except (NETWORK_ERRORS + (IOError,)) as e:
            log(f"  FAILED: {e}")
            failed += 1

    update_provenance(PROVENANCE, ACCESSION, render_table(list(manifest.values())))
    log(f"provenance written: {PROVENANCE}")
    log(f"summary: {ok} downloaded/recorded, {skipped} skipped, {failed} failed")

    if failed:
        log("NON-ZERO EXIT: some files failed. data/raw is unchanged for those.")
        return 1
    log("00_download.py complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
