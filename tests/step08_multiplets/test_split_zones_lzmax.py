"""Step 3 / LSRM Lzmax — unit для _split_zones_lzmax (3 кейса).

Цензор 2026-06-21: коммит 0c8726f заявлял «синтетический unit, 5 кейсов,
all pass», но теста в suite не было (false-green). Этот файл закрывает
блокер: (a) over-long zone дробится в долине min-counts, (b) зеркало
~7.47 ПШПВ (как реальный M3 Ac-228 233/252/277) НЕ дробится,
(c) исчерпание `_max_depth` наблюдаемо через `RuntimeWarning`.

Источник правила: pdf.md:482-485, 490-494 (LSRM Гамма-1С Lzmax).
Архитектурное обоснование длины зоны и порога: CLAUDE.md, раздел
«Ширина зоны / зонирование спектра».
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from gamma.peaks.deconvolve import _split_zones_lzmax  # noqa: E402


class _FakeLM:
    """Минимальный stand-in для LineMatch (только peak_channel)."""

    def __init__(self, ch: float, nuclide: str = "X", E: float = 0.0):
        self.peak_channel = float(ch)
        self.nuclide = nuclide
        self.library_E_keV = float(E)


def _flat_fwhm(value: float):
    return lambda ch: float(value)


def _make_spec(n_ch: int, valley_at: int, base: float = 1000.0,
               valley: float = 50.0):
    """Спектр с одной выраженной долиной в `valley_at`."""
    counts = np.full(n_ch, base, dtype=np.float64)
    counts[valley_at] = valley
    return SimpleNamespace(counts=counts)


# ──────────────────────────────────────────────────────────────────
# (a) over-long zone splits in min-counts valley
# ──────────────────────────────────────────────────────────────────

def test_a_overlong_zone_splits_at_min_counts_valley():
    """ПШПВ=10, пики на каналах 100 и 250 → размах=150, wing=2.5·10=25,
    zone_len = 150+50 = 200 > 10·10 = 100 → должна делиться.
    Долина min-counts задана на канале 175 → split-канал = 175.
    Ожидание: 2 группы, левая до 175 включительно, правая после.
    """
    members = [_FakeLM(100.0, "A", 100.0), _FakeLM(250.0, "B", 250.0)]
    clusters = [members]
    spec = _make_spec(n_ch=400, valley_at=175)

    out = _split_zones_lzmax(
        clusters, spec, _flat_fwhm(10.0),
        max_zone_length_fwhm=10.0, roi_window_factor=2.5,
    )

    assert len(out) == 2, (
        f"over-long zone должна была разбиться на 2 подзоны; "
        f"получено {len(out)}: {[[m.peak_channel for m in g] for g in out]}"
    )
    left, right = out[0], out[1]
    assert [m.peak_channel for m in left] == [100.0], (
        f"левая подзона должна содержать только канал 100 (split на 175); "
        f"got {[m.peak_channel for m in left]}"
    )
    assert [m.peak_channel for m in right] == [250.0], (
        f"правая подзона должна содержать только канал 250 (split на 175); "
        f"got {[m.peak_channel for m in right]}"
    )


# ──────────────────────────────────────────────────────────────────
# (b) ~7.47 FWHM зона (как реальный M3 Ac-228 233/252/277) НЕ дробится
# ──────────────────────────────────────────────────────────────────

def test_b_zone_about_747_fwhm_not_split():
    """ПШПВ=10, пики на каналах 100, 122, 150 → размах=50, wing=2·25=50,
    zone_len = 50+50 = 100. Порог 10·ПШПВ = 100 → '>' СТРОГИЙ → НЕ делится.
    Дополнительно даём долину между крайними пиками — split-путь не должен
    активироваться вообще (граница порога).
    """
    members = [
        _FakeLM(100.0, "Ac-228", 233.0),
        _FakeLM(122.0, "Ac-228", 252.0),
        _FakeLM(150.0, "Ac-228", 277.0),
    ]
    clusters = [members]
    spec = _make_spec(n_ch=300, valley_at=130)

    out = _split_zones_lzmax(
        clusters, spec, _flat_fwhm(10.0),
        max_zone_length_fwhm=10.0, roi_window_factor=2.5,
    )

    assert len(out) == 1, (
        f"зона ≈7.5·ПШПВ (как M3) НЕ должна дробиться при пороге 10·ПШПВ; "
        f"got {len(out)} подзон: "
        f"{[[m.peak_channel for m in g] for g in out]}"
    )
    assert [m.peak_channel for m in out[0]] == [100.0, 122.0, 150.0], (
        f"состав зоны должен совпадать со входом; "
        f"got {[m.peak_channel for m in out[0]]}"
    )


# ──────────────────────────────────────────────────────────────────
# (c) _max_depth exhaustion observability (RuntimeWarning)
# ──────────────────────────────────────────────────────────────────

def test_c_max_depth_exhaustion_emits_runtime_warning():
    """Та же over-long зона, что в (a), но `_max_depth=0` → split-путь НЕ
    исполняется (исчерпан лимит). Censor требует: это не должно молча
    проглатываться. Ожидание: ровно один RuntimeWarning с подстрокой
    `_max_depth`, и зона возвращается единым куском.
    """
    members = [_FakeLM(100.0, "A", 100.0), _FakeLM(250.0, "B", 250.0)]
    clusters = [members]
    spec = _make_spec(n_ch=400, valley_at=175)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = _split_zones_lzmax(
            clusters, spec, _flat_fwhm(10.0),
            max_zone_length_fwhm=10.0, roi_window_factor=2.5,
            _max_depth=0,
        )

    rw = [w for w in caught if issubclass(w.category, RuntimeWarning)
          and "_max_depth" in str(w.message)]
    assert len(rw) == 1, (
        f"ожидался ровно один RuntimeWarning об исчерпании _max_depth; "
        f"got {len(rw)} (всего пойманных: {len(caught)}): "
        f"{[str(w.message) for w in caught]}"
    )
    assert len(out) == 1 and len(out[0]) == 2, (
        f"при исчерпанном _max_depth зона должна вернуться единым куском; "
        f"got {[[m.peak_channel for m in g] for g in out]}"
    )