"""Verify every file in data/raw/ against its recorded SHA256."""
import hashlib, json, pathlib
ok = bad = 0
for mf in sorted(pathlib.Path('data/raw').rglob('.download_manifest.json')):
    for name, rec in json.loads(mf.read_text()).items():
        f = mf.parent / name
        if not f.exists():
            print(f"  MISSING  {f}"); bad += 1; continue
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        if h == rec['sha256']:
            ok += 1
        else:
            print(f"  MISMATCH {f}\n    recorded {rec['sha256']}\n    on disk  {h}")
            bad += 1
nf = pathlib.Path('data/raw/neftel_signatures/NIHMS1532254-supplement-9.xlsx')
exp = "208e73ab3d22c494caf85c867d69dc6be38df3fc62ab1f043d7fcc5441066277"
if nf.exists():
    h = hashlib.sha256(nf.read_bytes()).hexdigest()
    good = h == exp
    print(f"  {'OK      ' if good else 'MISMATCH'} {nf.name} (manual input)")
    ok += good; bad += (not good)
else:
    print("  MISSING  neftel manual input"); bad += 1
print(f"  raw files verified: {ok} ok, {bad} bad")
raise SystemExit(1 if bad else 0)
