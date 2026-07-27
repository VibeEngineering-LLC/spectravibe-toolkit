"""
Ollama-delegated generation: pytest suite for the 4 new diagnostic fields in
LITE BackgroundSubtractionResult (bg_subtract_energy.py, Phase 2 output).

POSTs prompt with current bg_subtract_energy.py source + _make_spec helper
precedent from the FULL safety test to qwen3-coder:30b at
http://127.0.0.1:11434/api/generate with format='json'. Saves the generated
pytest file into _drafts/bg_phase3/test_bg_subtract_energy_v1.py for human
review before promotion to tests/snapshot/.

Idempotent: re-running overwrites the draft. No production file is modified.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "scripts" / "gamma" / "calibration" / "bg_subtract_energy.py"
DRAFT_DIR = REPO / "_drafts" / "bg_phase3"
DRAFT = DRAFT_DIR / "test_bg_subtract_energy_v1.py"
RAW_LOG = DRAFT_DIR / "ollama_raw_response.json"

MODEL = "qwen3-coder:30b"
URL = "http://127.0.0.1:11434/api/generate"

SYSTEM_PROMPT = """You are a precise Python test-author engine.
You receive a Python module source and a SPEC describing pytest tests to generate.
You output STRICT JSON with one field: "file_content" containing the full pytest file.
Rules:
- Output is a single self-contained pytest file (Python 3.10+).
- Use numpy as np and import pytest.
- Reuse the EXACT _make_spec helper signature from the precedent block.
- Each test is a top-level function starting with `test_`.
- Each test has a short docstring stating WHAT is verified.
- Use pytest.approx for float comparisons, rel=1e-9 for exact-math fields.
- Use math.isnan or np.isnan for NaN-checks, never `== nan`.
- Do NOT add fixtures unless strictly required; per-test setup inline is fine.
- Do NOT add markers, parametrize, or conftest imports.
- Return ONLY the JSON object. No markdown fences, no commentary.
"""

SPEC = """SPEC for tests/snapshot/test_bg_subtract_energy.py.

Module under test: gamma.calibration.bg_subtract_energy
Public API:
    from gamma.calibration.bg_subtract_energy import (
        subtract_background, BackgroundSubtractionResult,
    )
Signature:
    subtract_background(sample, background, *, clamp_negative_to_zero=True) -> BackgroundSubtractionResult
Note: LITE variant — there is NO `user_confirmed_applicable` kwarg (that is FULL only).

Helper (copy verbatim, do not modify):

    from gamma.spectrum import Spectrum

    def _make_spec(counts, live_time, *, a0=0.0, a1=1.0, source_path="test"):
        spec = Spectrum(
            counts=np.array(counts, dtype=np.int64),
            live_time=float(live_time),
            real_time=float(live_time),
            source_path=source_path,
            source_format="test",
        )
        spec.energy_cal = (float(a0), float(a1))
        spec.n_channels = len(counts)
        spec.n_channels_raw = len(counts)
        return spec

Generate exactly NINE pytest functions in this order:

1. test_n_channels_clipped_counts_negative_diff
   - sample = [10, 10, 10, 10, 10] live_time=100, a0=0, a1=1
   - bg     = [ 5,  5, 20, 25, 30] live_time=100, a0=0, a1=1
   - Call subtract_background(sample, bg) (default clamp=True).
   - With k=1.0, raw_diff = [5, 5, -10, -15, -20] -> 3 channels negative.
   - Assert result.n_channels_clipped == 3.

2. test_n_channels_clipped_zero_when_sample_dominates
   - sample = [100]*5 live_time=100, a0=0, a1=1
   - bg     = [1]*5   live_time=100, a0=0, a1=1
   - Assert result.n_channels_clipped == 0.

3. test_net_uncertainties_poisson_propagation_equal_live_times
   - sample = [100]*4 live_time=100, a0=0, a1=1
   - bg     = [25]*4  live_time=100, a0=0, a1=1
   - k=1.0; expected sigma_net = sqrt(100 + 1*25) = sqrt(125) per channel.
   - Assert result.net_uncertainties is not None.
   - Assert result.net_uncertainties.shape == (4,)
   - Assert each element == pytest.approx(np.sqrt(125.0), rel=1e-9).

4. test_net_uncertainties_with_scale_factor_k_neq_1
   - sample = [200]*4 live_time=200, a0=0, a1=1
   - bg     = [50]*4  live_time=100, a0=0, a1=1
   - k = 200/100 = 2.0; bg_on_sample = 50*2 = 100; expected sigma = sqrt(200 + 2*100) = sqrt(400) = 20.0
   - Assert result.scale_factor == pytest.approx(2.0, rel=1e-12).
   - Assert each net_uncertainties element == pytest.approx(20.0, rel=1e-9).

5. test_gain_mismatch_relative_known_a1
   - sample = [10]*8 live_time=100, a0=0.0, a1=3.0
   - bg     = [10]*8 live_time=100, a0=0.0, a1=2.97
   - Expected gain_mismatch_relative = abs(3.0 - 2.97) / max(3.0, 2.97) = 0.03/3.0 = 0.01
   - Assert result.gain_mismatch_relative == pytest.approx(0.01, rel=1e-9).

6. test_gain_mismatch_relative_nan_when_a1_missing
   - Build a normal sample and bg with a1=1.0 first via _make_spec.
   - Then MUTATE one of them: bg.energy_cal = (0.0,) (tuple of length 1, no a1 index).
   - subtract_background should still run; result.gain_mismatch_relative must be NaN.
   - Use math.isnan or np.isnan for the assertion (import math at top).

7. test_zero_point_mismatch_keV_known_a0
   - sample = [10]*8 live_time=100, a0=0.0, a1=3.0
   - bg     = [10]*8 live_time=100, a0=4.5, a1=3.0
   - Expected zero_point_mismatch_keV = abs(0.0 - 4.5) = 4.5
   - Assert result.zero_point_mismatch_keV == pytest.approx(4.5, rel=1e-9).
   - Also assert F-243 note phrase is NOT in result.notes (4.5 <= 5.0 threshold).

8. test_notes_extension_above_5_keV_threshold
   - sample = [10]*8 live_time=100, a0=0.0, a1=3.0
   - bg     = [10]*8 live_time=100, a0=10.0, a1=3.0  (Delta a0 = 10 keV > 5)
   - Assert result.zero_point_mismatch_keV == pytest.approx(10.0, rel=1e-9).
   - Assert the substring "F-243" is in result.notes.
   - Assert the substring "bg_subtract_dual_mode.py" is in result.notes.
   - Assert the substring "keV" is in result.notes.

9. test_clamp_negative_false_preserves_negative_net_but_count_unchanged
   - sample = [10, 10, 10, 10] live_time=100, a0=0, a1=1
   - bg     = [ 1, 50, 50,  1] live_time=100, a0=0, a1=1
   - Call subtract_background(sample, bg, clamp_negative_to_zero=False).
   - Expected raw diff = [9, -40, -40, 9].
   - Assert result.n_channels_clipped == 2.
   - Assert result.net_counts contains at least one negative value (any(result.net_counts < 0)).
   - Then call again with clamp_negative_to_zero=True (default): assert min(result.net_counts) >= 0.0
     AND n_channels_clipped is still 2 (counted BEFORE clamp).

GENERAL:
- File header docstring: one short paragraph stating purpose ("Snapshot tests for the 4 diagnostic fields added to BackgroundSubtractionResult in Phase 2 of the F-58 LITE/FULL untangle.").
- Imports at top: from __future__ import annotations; import math; import numpy as np; import pytest; from gamma.spectrum import Spectrum; from gamma.calibration.bg_subtract_energy import subtract_background, BackgroundSubtractionResult.
- _make_spec helper directly under imports (before the test functions).
- No extra helpers, no parametrize, no class wrappers.
- Return ONLY {"file_content": "<full text>"}.
"""


def call_ollama(current_src: str) -> dict:
    user_prompt = (
        f"SPEC:\n{SPEC}\n\n"
        f"MODULE UNDER TEST (verbatim, for context only, do not echo):\n```python\n{current_src}\n```\n\n"
        "Return JSON: {\"file_content\": \"<full new pytest file text>\"}\n"
    )
    payload = {
        "model": MODEL,
        "prompt": SYSTEM_PROMPT + "\n\n" + user_prompt,
        "stream": False,
        "format": "json",
        "options": {"temperature": 0.0, "num_predict": 8192},
    }
    r = requests.post(URL, json=payload, timeout=600)
    r.raise_for_status()
    return r.json()


def main() -> int:
    DRAFT_DIR.mkdir(parents=True, exist_ok=True)
    current = SRC.read_text(encoding="utf-8")
    print(f"[info] module under test {SRC.name}: {len(current)} chars, {current.count(chr(10))+1} lines", flush=True)

    print(f"[info] POST {URL} model={MODEL}, this can take 60-300s ...", flush=True)
    resp = call_ollama(current)
    RAW_LOG.write_text(json.dumps(resp, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[info] raw ollama response saved -> {RAW_LOG}", flush=True)

    inner = resp.get("response", "")
    if not inner:
        print("[ERR] empty response field; see raw log", file=sys.stderr)
        return 2
    try:
        parsed = json.loads(inner)
    except json.JSONDecodeError as e:
        print(f"[ERR] response is not JSON: {e}; first 500 chars: {inner[:500]!r}", file=sys.stderr)
        return 3
    fc = parsed.get("file_content")
    if not isinstance(fc, str) or len(fc) < 200:
        print(f"[ERR] file_content missing or suspiciously short: {len(fc) if isinstance(fc, str) else type(fc).__name__}", file=sys.stderr)
        return 4

    DRAFT.write_text(fc, encoding="utf-8")
    print(f"[OK] draft saved -> {DRAFT} ({len(fc)} chars, {fc.count(chr(10))+1} lines)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())