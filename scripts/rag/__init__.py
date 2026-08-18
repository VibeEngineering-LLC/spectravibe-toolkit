# -*- coding: utf-8 -*-
"""scripts/rag — RAG index generators for SpectraVibe.

Houses:
- `build_spectra_index.py` — mechanically-derived spectrum-RAG
  (`audit/_rag/SPECTRA_INDEX.json`), regenerated from an operator's
  LSRM `.spe` tree. Separate from methodology-RAG (`RAG_INDEX.json`,
  hand-curated).

F-150 / F-115 compliance: paths in the emitted index use `<LSRM>`
placeholder; no operator-absolute paths leaked into the repo.
"""
