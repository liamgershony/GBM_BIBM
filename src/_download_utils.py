"""Shared download helpers for 00_download.py and 00b_download_cgga.py.

Standard library only, on purpose: these scripts must run before
envs/environment.yml exists and before any conda environment is built.

Nothing here parses, decompresses or loads scientific data. Files are fetched as
bytes, hashed, and recorded. data/raw/ is immutable (CLAUDE.md 7.4) so a file
that already exists and matches its recorded checksum is never rewritten.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

USER_AGENT = "gbm-persister/0.1 (academic use; contact liam.gershony@gmail.com)"
CHUNK = 1024 * 1024  # 1 MiB

# The only failures a network call is allowed to swallow. Anything else -- a
# NameError, an AttributeError, a ModuleNotFoundError -- is a bug in our code and
# must propagate. A broad `except Exception` here turns a programming error into a
# silent "server unavailable" and the pipeline reports success.
NETWORK_ERRORS = (urllib.error.URLError, urllib.error.HTTPError,
                  TimeoutError, OSError, ConnectionError)


def _ssl_context() -> ssl.SSLContext:
    """A fully-verifying TLS context that works on python.org macOS installs.

    The python.org installer ships its own trust store and leaves it empty until
    you run "Install Certificates.command", so urllib fails with
    CERTIFICATE_VERIFY_FAILED on sites curl handles fine. We look for a real CA
    bundle in preference order and use it.

    Certificate verification is NEVER disabled. These files are the provenance
    record for the whole paper; fetching them over an unverified channel would
    make every SHA256 below meaningless.
    """
    try:
        import certifi  # optional; present once envs/environment.yml is built
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        pass
    for bundle in ("/etc/ssl/cert.pem", "/usr/local/etc/openssl/cert.pem"):
        if os.path.exists(bundle):
            return ssl.create_default_context(cafile=bundle)
    return ssl.create_default_context()


SSL_CONTEXT = _ssl_context()


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

def make_logger(logfile: Path):
    """Return log(msg). Writes to stdout and appends to logfile, line-buffered,
    so a background run can be followed with `tail -f`."""
    logfile.parent.mkdir(parents=True, exist_ok=True)
    fh = open(logfile, "a", buffering=1, encoding="utf-8")

    def log(msg: str) -> None:
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        line = f"[{stamp}] {msg}"
        print(line, flush=True)
        fh.write(line + "\n")

    return log


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------
# http
# --------------------------------------------------------------------------

def _request(url: str, headers: dict | None = None):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return req


def fetch_text(url: str, timeout: int = 60, retries: int = 4, log=print) -> str:
    """GET a small text/HTML resource (a directory index or download page)."""
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(_request(url), timeout=timeout, context=SSL_CONTEXT) as r:
                return r.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            if attempt == retries:
                raise
            wait = 2 ** attempt
            log(f"  listing fetch failed ({e}); retry {attempt}/{retries - 1} in {wait}s")
            time.sleep(wait)
    raise RuntimeError("unreachable")


def remote_size(url: str, timeout: int = 60) -> int | None:
    """Content-Length via HEAD, or None if the server does not report it."""
    req = _request(url)
    req.get_method = lambda: "HEAD"
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as r:
            n = r.headers.get("Content-Length")
            if n is None:
                return None
            # A HEAD that answers Content-Length: 0 is telling us nothing useful --
            # some endpoints (Europe PMC's supplementaryFiles among them) do this
            # for generated payloads. Treating 0 as a real size makes an absent
            # .part look "already complete". Report unknown instead.
            return int(n) or None
    except NETWORK_ERRORS:
        # HEAD is advisory: callers treat None as "size unknown" and skip size
        # verification. Only network failures may produce that; anything else
        # would silently disable a correctness check.
        return None
    except ValueError:
        return None                      # non-integer Content-Length


def extract_links(html: str, base_url: str) -> list[tuple[str, str]]:
    """Return (filename, absolute_url) for every href in an index page.

    Skips parent-directory links, sort links, and sub-directories. Deliberately
    dumb: we do not guess filenames, we read whatever the server lists.
    """
    out: list[tuple[str, str]] = []
    seen: set[str] = set()
    for href in re.findall(r'href\s*=\s*["\']([^"\']+)["\']', html, flags=re.I):
        if href.startswith(("?", "#", "mailto:")) or href in ("../", "/"):
            continue
        if href.endswith("/"):
            continue
        name = href.rsplit("/", 1)[-1]
        if not name or name in seen:
            continue
        seen.add(name)
        url = href if href.startswith(("http://", "https://")) else base_url.rstrip("/") + "/" + href.lstrip("/")
        out.append((name, url))
    return out


# --------------------------------------------------------------------------
# identifier normalisation
# --------------------------------------------------------------------------

# Tokens that mean "no identifier" when they appear in an ID or grouping column.
# GSE174554's GEO `pair#` and Wang et al. Supplementary Table 1's `Pair#` BOTH use
# the literal string "NA" for unpaired specimens.
NULL_ID_TOKENS = frozenset({"", "NA", "N/A", "NONE", "NULL", "-", "."})


def read_id_column(value) -> str | None:
    """Normalise one identifier cell to a real id, or None if it means 'absent'.

    Use this for EVERY column whose values are used to group or join records --
    patient ids, pair numbers, specimen ids. Never compare a raw cell directly.

    Why this exists. A literal "NA" left in a grouping column is not inert: it is a
    valid dict key, so every unpaired specimen collapses into one group. During
    cohort exploration this fabricated a patient TWICE, in two separate files:

      * GEO `pair#`      -- 15 specimens fused into one "patient" that appeared to
                            hold 3 primaries and 9 recurrents, inflating the
                            matched-pair count from 30 to 31.
      * Supplementary `Pair#` -- 10 specimens fused the same way, presenting as a
                            pair the GEO scheme had "lost" and again giving 31.

    Both times the corruption surfaced as a plausible *finding* rather than as an
    error, which is exactly why it must be handled in one shared place instead of
    at each call site.

    Returns None (not "") so that a missing id is falsy AND distinguishable from a
    genuine empty-string id, and so csv writers emit a blank rather than "None".
    """
    if value is None:
        return None
    s = str(value).strip()
    return None if s.upper() in NULL_ID_TOKENS else s


# --------------------------------------------------------------------------
# hashing
# --------------------------------------------------------------------------

# Magic bytes by extension. A server that answers 200 with an HTML interstitial
# instead of the file (PMC's anti-bot proof-of-work page does exactly this) would
# otherwise be saved under the requested filename and checksummed, turning a
# failed download into a recorded "success". This is the same failure shape as the
# Content-Length: 0 bug and the batch2 false joins.
MAGIC = {
    ".xlsx": (b"PK",), ".zip": (b"PK",), ".gz": (b"\x1f\x8b",),
    ".tar": (b"", ), ".pdf": (b"%PDF",),
}


def verify_file_type(path: Path) -> None:
    """Raise if the bytes on disk do not match what the extension promises."""
    exts = [e for e in MAGIC if path.name.lower().endswith(e)]
    if not exts:
        return
    expect = MAGIC[exts[0]]
    if not any(expect):
        return
    with open(path, "rb") as fh:
        head = fh.read(8)
    if not any(head.startswith(m) for m in expect if m):
        looks_html = b"<html" in head.lower() or head.startswith(b"<!DO")
        raise IOError(
            f"{path.name}: content does not match its type -- expected one of "
            f"{[m for m in expect]}, got {head!r}"
            + (". The server returned an HTML page (login wall, anti-bot "
               "challenge, or error) instead of the file." if looks_html else ""))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


# --------------------------------------------------------------------------
# resumable download
# --------------------------------------------------------------------------

def download_with_resume(url: str, dest: Path, timeout: int = 120,
                         retries: int = 5, log=print) -> int:
    """Download url -> dest, resuming a partial .part file via HTTP Range.

    Returns the final byte count. Writes to dest.part and renames atomically
    only after the transfer completes, so dest never exists in a partial state.
    """
    part = dest.with_suffix(dest.suffix + ".part")
    expected = remote_size(url, timeout=timeout)

    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0

        if expected is not None and have == expected:
            log(f"  .part already complete ({have:,} B)")
            break
        if expected is not None and have > expected:
            log(f"  .part larger than remote ({have:,} > {expected:,}) — restarting")
            part.unlink()
            have = 0

        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with urllib.request.urlopen(_request(url, headers), timeout=timeout, context=SSL_CONTEXT) as r:
                # 206 = server honoured Range; 200 = it ignored it, so start over.
                if have and r.status == 200:
                    log("  server ignored Range; restarting from 0")
                    have = 0
                mode = "ab" if have else "wb"
                with open(part, mode) as f:
                    while True:
                        block = r.read(CHUNK)
                        if not block:
                            break
                        f.write(block)
                    f.flush()
                    os.fsync(f.fileno())
            got = part.stat().st_size
            if expected is None or got == expected:
                break
            log(f"  short read ({got:,} of {expected:,} B) — resuming")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            if attempt == retries:
                raise
            wait = min(60, 2 ** attempt)
            log(f"  transfer error ({e}); retry {attempt}/{retries - 1} in {wait}s")
            time.sleep(wait)

    if not part.exists():
        raise IOError(f"download produced no data for {dest.name} "
                      f"(no {part.name} on disk after {retries} attempt(s))")
    final = part.stat().st_size
    if expected is not None and final != expected:
        # Keep the .part for inspection. Do not record a checksum for a file we
        # cannot vouch for.
        raise IOError(f"size mismatch for {dest.name}: got {final:,} B, expected {expected:,} B "
                      f"(partial kept at {part})")
    part.replace(dest)
    try:
        verify_file_type(dest)
    except IOError:
        # Keep the evidence but never leave a mislabelled file in data/raw/.
        bad = dest.with_suffix(dest.suffix + ".rejected")
        dest.replace(bad)
        raise
    return final


# --------------------------------------------------------------------------
# manifest  (our own record; GEO/CGGA do not publish checksums)
# --------------------------------------------------------------------------

def load_manifest(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_manifest(path: Path, manifest: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(path)


# --------------------------------------------------------------------------
# PROVENANCE.md
# --------------------------------------------------------------------------

def render_table(entries: list[dict]) -> str:
    if not entries:
        return "_No files recorded yet._"
    rows = ["| File | Accession | URL | Downloaded (UTC) | Bytes | SHA256 |",
            "|---|---|---|---|---|---|"]
    for e in sorted(entries, key=lambda x: x["file"]):
        rows.append(f"| `{e['file']}` | {e['accession']} | <{e['url']}> | "
                    f"{e['downloaded_utc']} | {e['bytes']:,} | `{e['sha256']}` |")
    return "\n".join(rows)


def update_provenance(path: Path, marker: str, block: str) -> None:
    """Replace the region between the BEGIN/END markers for one accession.

    Idempotent: rerunning a download rewrites only its own block and leaves
    every other section of PROVENANCE.md untouched.
    """
    begin = f"<!-- BEGIN AUTOGENERATED: {marker} -->"
    end = f"<!-- END AUTOGENERATED: {marker} -->"
    text = path.read_text()
    new_region = f"{begin}\n{block}\n{end}"
    if begin in text and end in text:
        pre = text.split(begin)[0]
        post = text.split(end, 1)[1]
        text = pre + new_region + post
    else:
        text = text.rstrip() + f"\n\n### {marker}\n\n{new_region}\n"
    path.write_text(text)
