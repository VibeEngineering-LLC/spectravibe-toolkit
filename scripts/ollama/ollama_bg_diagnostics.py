"""
Ollama-delegated upgrade: add 4 diagnostic fields to LITE BackgroundSubtractionResult.

POSTs prompt with current bg_subtract_energy.py to qwen3-coder:30b at
http://127.0.0.1:11434/api/generate with format='json'. Saves the regenerated file
into _drafts/bg_phase2/bg_subtract_energy_v2.py for human review.

Idempotent: re-running overwrites the draft. No production file is modified here.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import requests

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "scripts" / "gamma" / "calibration" / "bg_subtract_energy.py"
DRAFT_DIR = REPO / "_drafts" / "bg_phase2"
DRAFT = DRAFT_DIR / "bg_subtract_energy_v2.py"
RAW_LOG = DRAFT_DIR / "ollama_raw_response.json"

MODEL = "qwen3-coder:30b"
URL = "http://127.0.0.1:11434/api/generate"

SYSTEM_PROMPT = """You are a precise Python refactor engine.
You receive a Python source file and a SPEC describing additive changes.
You output STRICT JSON with one field: "file_content" containing the full new file.
Rules:
- Preserve every existing public name (dataclass, fields, function names, __all__).
- Preserve existing field order; APPEND new fields at the end.
- Do NOT add new imports beyond numpy / dataclasses / typing / __future__ unless strictly required.
- Do NOT change docstring of the module's first three paragraphs; you MAY add a paragraph at the end about new diagnostics.
- Do NOT alter the public signature of subtract_background.
- Do NOT add type-stub files, tests, or examples.
- Return ONLY the JSON object. No markdown fences, no commentary.
"""

SPEC = """SPEC вЂ” additive diagnostics for BackgroundSubtractionResult (file: bg_subtract_energy.py).

ADD exactly 4 new fields to the dataclass, AFTER existing fields, with these EXACT names, types, and defaults:

    net_uncertainties: np.ndarray | None = None        # sqrt(N_src + k^2 * N_bg) per channel, length = sample.n_channels
    gain_mismatch_relative: float = float("nan")       # |a1_sample - a1_bg| / max(|a1_sample|, |a1_bg|); nan if either energy_cal lacks index 1
    zero_point_mismatch_keV: float = float("nan")      # abs(a0_sample - a0_bg); nan if either energy_cal lacks index 0
    n_channels_clipped: int = 0                        # count of channels where raw_diff = sample_counts - bg_on_sample < 0 (BEFORE clamp)

COMPUTATION inside subtract_background():
- n_channels_clipped: count of channels where (sample_counts - bg_on_sample) < 0 (use np.sum on a boolean mask, BEFORE any clamping).
- net_uncertainties: per-channel Poisson sigma propagation: sigma_net[i] = sqrt(sample_counts[i] + (k**2) * bg_on_sample_unscaled[i]).
  NOTE: bg_on_sample is ALREADY scaled by k in the existing code (the multiplication * k happens inside the np.interp expression). For the uncertainty you need the UNSCALED interpolated bg-counts; compute them as bg_on_sample / k (or recompute the interp separately as bg_on_sample_raw before scaling). Be careful: if k > 0 always (asserted above), division is safe; use np.where(k > 0, bg_on_sample / k, 0.0) for defensiveness. Equivalent identity: sqrt(sample_counts + k * bg_on_sample). Either form is acceptable; prefer the second (sample_counts + k * bg_on_sample) because it avoids the division step entirely. Mathematical equivalence: k^2 * (bg/k) = k * bg, so sqrt(N_src + k^2 * (bg/k)) = sqrt(N_src + k * bg). Use this simpler form.
- gain_mismatch_relative: extract a1_s = sample.energy_cal[1] if len(sample.energy_cal) >= 2 else None; same for a1_b. If either is None or max(|a1_s|, |a1_b|) == 0 => float('nan'). Else abs(a1_s - a1_b) / max(abs(a1_s), abs(a1_b)).
- zero_point_mismatch_keV: extract a0_s = sample.energy_cal[0] if len(sample.energy_cal) >= 1 else None; same for a0_b. If either is None => float('nan'). Else abs(a0_s - a0_b).

NOTES STRING extension:
- After the existing overlap warning, append (only if zero_point_mismatch_keV is a finite number AND > 5.0):
  " О”aв‚Ђ={zero_point_mismatch_keV:.2f} keV (СЃРј. F-243 РІ bg_subtract_dual_mode.py РґР»СЏ strict-mode subtraction)."
- The О”aв‚Ђ phrase MUST use the literal Cyrillic О” and lowercase greek 0. The keV unit MUST be ASCII 'keV'.

FIELDS in BackgroundSubtractionResult constructor call at the end of the function MUST be passed by name. Pre-existing fields keep their current values; new fields are populated from local variables you compute above.

Do not change __all__.
Do not introduce a separate helper function unless absolutely needed.
"""


def call_ollama(current_src: str) -> dict:
    user_prompt = (
        f"SPEC:\n{SPEC}\n\nCURRENT FILE CONTENT (verbatim):\n```python\n{current_src}\n```\n\n"
        "Return JSON: {\"file_content\": \"<full new file text>\"}\n"
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
    print(f"[info] current file {SRC.name}: {len(current)} chars, {current.count(chr(10))+1} lines", flush=True)

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
