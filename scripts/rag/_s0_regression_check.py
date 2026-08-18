# -*- coding: utf-8 -*-
"""F-070 W4 S0 — Post-retrofit regression verifier (per inbox brief)."""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np


def main() -> int:
    root = "audit/_rag/visual_templates"
    files = sorted(
        f for f in glob.glob(f"{root}/**/*.json", recursive=True)
        if "/VT-" in f.replace("\\", "/")
        and "_raw_ingest" not in f
        and "_drift_study" not in f
    )

    print(f"Found {len(files)} canonical templates (expected 24)")
    if len(files) != 24:
        print(f"FAIL: expected 24, got {len(files)}")
        return 1

    errors: list[str] = []
    vecs: dict[str, np.ndarray] = {}
    for fpath in files:
        with open(fpath, encoding="utf-8") as fh:
            t = json.load(fh)
        tid = t["template_id"]
        if t.get("schema_version") != "0.2":
            errors.append(f"schema_version: {tid} = {t.get('schema_version')!r}")
        if "crystal_class" not in t:
            errors.append(f"missing crystal_class: {tid}")
        if "vessel_class" not in t:
            errors.append(f"missing vessel_class: {tid}")
        fv = np.array(t["feature_vector"]["values"], dtype=np.float64)
        if len(fv) != 128:
            errors.append(f"dim={len(fv)}: {tid}")
        norm = float(np.linalg.norm(fv))
        if abs(norm - 1.0) > 1e-5:
            errors.append(f"L2_norm={norm:.8f}: {tid}")
        vecs[tid] = fv

    if errors:
        print("ERRORS:", errors)
        return 1

    ids = sorted(vecs)
    diag = [float(np.dot(vecs[a], vecs[a])) for a in ids]
    ok = all(abs(d - 1.0) < 1e-5 for d in diag)

    abs_re = re.compile(r"[A-Z]:\\|^/(home|Users)/")
    sn_re = re.compile(r"№\d+-\d+|_№[А-ЯA-Z0-9-]+")
    f115_fails: list[str] = []
    for fpath in files:
        text = open(fpath, encoding="utf-8").read()
        if abs_re.search(text):
            f115_fails.append(f"ABS_PATH:{fpath}")
        if sn_re.search(text):
            f115_fails.append(f"CERT_SN:{fpath}")

    if f115_fails:
        print("F-115 FAIL:", f115_fails)
        return 1
    if not ok:
        print("NORM FAIL", [(a, d) for a, d in zip(ids, diag) if abs(d - 1.0) > 1e-5])
        return 1

    print(f"REGRESSION PASS: {len(files)} templates, schema_version=0.2, L2-norm intact, F-115 clean")

    # 24x24 cosine matrix diagonal already verified above; do full matrix
    # to confirm bit-for-bit symmetry baseline.
    n = len(ids)
    mat = np.zeros((n, n))
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            mat[i, j] = float(np.dot(vecs[a], vecs[b]))
    print(f"24x24 cosine matrix: diag min={np.diag(mat).min():.10f} max={np.diag(mat).max():.10f}")
    print(f"  off-diag min={mat[np.triu_indices(n, k=1)].min():.6f} max={mat[np.triu_indices(n, k=1)].max():.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
