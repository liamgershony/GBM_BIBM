# PROVENANCE — data/raw/

`data/raw/` is **immutable**. Never write to it, never edit in place. Every
downstream file must be regenerable from `data/raw/` + `src/` + `configs/`.
(CLAUDE.md §7.4)

One row per downloaded file. Fill in at download time, not afterwards.

| File | Source / accession | URL | Download date (UTC) | Bytes | SHA256 |
|---|---|---|---|---|---|
| _(pending)_ | GSE174554 | | | | |
| _(pending)_ | CGGA mRNAseq_693 | | | | |
| _(pending)_ | CGGA mRNAseq_325 | | | | |

All three sources are open access: no account, no data-use agreement, no approval
queue. Verified 2026-08-23.

**Not used:** GLASS (Synapse account + DUA), GSE174554 exome (not in the public
deposit), spatial subcohort (out of scope).
