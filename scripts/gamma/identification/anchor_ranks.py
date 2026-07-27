"""
Practical anchor-peak ranking on NaI/CsI/RadiaCode (F-79 / v1.11.1).

User methodology (15.11.2025):
"Сортировка идёт от наиболее легко узнаваемых к наиболее неоднозначным."

The standard `identify_nuclides` characteristic-line picker is
"max-intensity line" — but on moderate-resolution detectors with
Compton-dominated background, the highest-intensity library line is
often NOT the most reliable for identification:
  • Pb-212 238.6 (I=43.6%) sits in a crowded 240-keV merged region
    with Pb-214 242 and Ra-224 241 on NaI 63×63 (FWHM ~24 keV).
  • Bi-214 609.3 (I=45.5%) sits in dense Compton continuum and
    competes with Cs-134 604 and other lines.
  • Tl-208 583.2 (I=30.6%) is often masked by Bi-214 609 on NaI.

The expert heuristic uses a different priority: anchors ranked by
*observed visibility* in real spectra under moderate resolution.

Rank | Energy   | Nuclide     | Note
-----|----------|-------------|-----------------------------------------
  1  | 2614.51  | Tl-208      | Эталонный маяк, почти нет конкурентов
  2  | 1460.82  | K-40        | Чистая область, мало конкурентов
  3a | 1173.23  | Co-60       | Парный паттерн (доказывается вместе с 3b)
  3b | 1332.49  | Co-60       | Парный паттерн (доказывается вместе с 3a)
  4  |  661.66  | Cs-137      | Яркий одиночный в чистой зоне
  5  |   59.54  | Am-241      | Характерный низкоэнергетический
  6  | 1764.49  | Bi-214      | Хорошо отделён от соседей
  7  |  911.20  | Ac-228      | Довольно чистая область
  8  |  609.31  | Bi-214      | Сильный, но в дозамещении Compton-полосой
  9  |  351.93  | Pb-214      | Узнаваем, зона перегружена
 10  |  295.22  | Pb-214      | Подтверждение радонового ряда
 11  |  583.19  | Tl-208      | Может теряться на фоне Bi-214
 12  |  238.63  | Pb-212      | Часто замаскирован Pb-XRF/Ra-224
 13  |  511.0   | annihilation| Неспецифичен — почти бесполезен в одиночку
 14  |  186.2   | U-235/Ra-226| Классическое перекрытие (HPGe-only resolution)

Express patterns (must be confirmed together):
  • Co-60        : 1173 + 1332 (ratio ≈ 1:1)
  • Bi-214       : 609 + 1764 (Ra chain confirmation pair)
  • Bi-214 quartet: 609 + 1120 + 1764 + 2204 (strong evidence)
  • Pb-214       : 352 + 295 (Pb-214 doublet)
  • Th-232 triplet: 2615 + 911 + 583 (Tl-208 / Ac-228 / Tl-208)
  • Th-232 strong : 2615 + 911 (minimum for Th)
  • Cs-137 pair  : 662 + Ba Kα 32.07 (IC daughter)

The detection algorithm:
  Pass A — find anchors ranks 1-7 directly by energy match
  Pass B — confirm express patterns (boosts CI of confirmed parents)
  Pass C — note "express identification" on the matched anchors
  Pass D — disputed regions below 400 keV go to identify_nuclides
            with full disambiguate Rules 1-5

This module provides the data tables and the matchers. The orchestrator
in `staged_pipeline.analyze_lsrm_spe` does the multi-pass logic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from gamma.data.nuclide_library import get_nuclide


# ──────────────────────────────────────────────────────────────────
# Anchor rank table (data)
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class AnchorEntry:
    rank: int
    energy_keV: float
    nuclide: str                    # "" for non-nuclide (annihilation)
    chain: str                      # "Th-232" / "U-238" / "" / "ambiguous"
    note: str
    requires_partner: bool = False  # True if must be confirmed with companion
    partner_energies: Tuple[float, ...] = ()


ANCHOR_RANKS: Tuple[AnchorEntry, ...] = (
    AnchorEntry(
        rank=1, energy_keV=2614.51, nuclide="Tl-208", chain="Th-232",
        note="Эталонный маяк, почти нет конкурентов"
    ),
    AnchorEntry(
        rank=2, energy_keV=1460.82, nuclide="K-40", chain="",
        note="Чистая область, мало природных конкурентов"
    ),
    AnchorEntry(
        rank=3, energy_keV=1173.23, nuclide="Co-60", chain="",
        note="Парная сигнатура с 1332.49 keV",
        requires_partner=True, partner_energies=(1332.49,)
    ),
    AnchorEntry(
        rank=3, energy_keV=1332.49, nuclide="Co-60", chain="",
        note="Парная сигнатура с 1173.23 keV",
        requires_partner=True, partner_energies=(1173.23,)
    ),
    AnchorEntry(
        rank=4, energy_keV=661.66, nuclide="Cs-137", chain="",
        note="Яркий одиночный пик в чистой зоне"
    ),
    AnchorEntry(
        rank=5, energy_keV=59.54, nuclide="Am-241", chain="",
        note="Характерный низкоэнергетический пик"
    ),
    AnchorEntry(
        rank=6, energy_keV=1764.49, nuclide="Bi-214", chain="U-238",
        note="Хорошо отделён от большинства линий (радоновый маркер)"
    ),
    AnchorEntry(
        rank=7, energy_keV=911.20, nuclide="Ac-228", chain="Th-232",
        note="Довольно чистая область — вторичный Th-маркер"
    ),
    AnchorEntry(
        rank=8, energy_keV=609.31, nuclide="Bi-214", chain="U-238",
        note="Сильный, но в области Compton-фона — подтверждать 1764"
    ),
    AnchorEntry(
        rank=9, energy_keV=351.93, nuclide="Pb-214", chain="U-238",
        note="Узнаваем, но зона перегружена"
    ),
    AnchorEntry(
        rank=10, energy_keV=295.22, nuclide="Pb-214", chain="U-238",
        note="Подтверждает радоновый ряд (Pb-214 doublet)"
    ),
    AnchorEntry(
        rank=11, energy_keV=583.19, nuclide="Tl-208", chain="Th-232",
        note="Может теряться на фоне Bi-214 609 на NaI"
    ),
    AnchorEntry(
        rank=12, energy_keV=238.63, nuclide="Pb-212", chain="Th-232",
        note="Часто слабый и замаскирован Pb-XRF/Pb-214"
    ),
    AnchorEntry(
        rank=13, energy_keV=511.0, nuclide="", chain="",
        note="Аннигиляция/Tl-208 510.77 — неспецифичен в одиночку"
    ),
    AnchorEntry(
        rank=14, energy_keV=186.2, nuclide="U-235/Ra-226", chain="ambiguous",
        note="Классическое перекрытие U-235 185.71 + Ra-226 186.21 — HPGe-only"
    ),
    # ── F-453-FU v2 calibration-tier anchors (2026-06-23) ──────────────
    # Закрывают BUG-38 high-E extrapolation на коротких NaI калибрационных
    # фикстурах (AmTiCsEu Marinelli). ДВУХУРОВНЕВЫЙ anti-shadow gate:
    #   (1) fixture-fingerprint gate — в `find_anchor_matches` calibration-
    #       tier пробуется ТОЛЬКО когда в peak list видны peak'и возле
    #       Cs-137 661.66 keV И Am-241 59.54 keV (подпись AmTiCsEu/AmCs);
    #       на Ra-226/Th-232/K-40/Co-60 — нет ни Cs ни Am, gate fail,
    #       calibration-tier silent.
    #   (2) partner-required (вторичная защита) — каждый anchor требует
    #       второго маркера фикстуры (Sc-44↔Ti-44 mutual; Eu-152 1408 →
    #       Sc-44 ИЛИ Ti-44, НЕ другие Eu-152 линии чтобы не получить
    #       false-positive через Ac-228 338/794).
    # Активируются через `find_anchor_matches` независимо от max_rank —
    # см. CALIBRATION_RANK_START и AMTICSEU_FINGERPRINT_LINES_KEV ниже.
    AnchorEntry(
        rank=15, energy_keV=67.87, nuclide="Sc-44", chain="",
        note="Sc-44 K-line (daughter Ti-44 equilibrium). Partner — Ti-44 1157",
        requires_partner=True, partner_energies=(1157.02,)
    ),
    AnchorEntry(
        rank=16, energy_keV=1157.02, nuclide="Ti-44", chain="",
        note="Ti-44 main γ. Partner — Sc-44 67.87 (daughter equilibrium)",
        requires_partner=True, partner_energies=(67.87,)
    ),
    AnchorEntry(
        rank=17, energy_keV=1408.01, nuclide="Eu-152", chain="",
        note=("Eu-152 high-E (BUG-38 high-E anchor). Partner — Sc-44 67.87 "
              "ИЛИ Ti-44 1157 (маркеры калибрационной фикстуры; не Eu-152 "
              "линии — те confound с Ac-228 338/794)."),
        requires_partner=True, partner_energies=(67.87, 1157.02)
    ),
)


# F-453-FU v2: ranks ≥ CALIBRATION_RANK_START — calibration-tier anchors,
# включаются в `find_anchor_matches` независимо от max_rank, но ТОЛЬКО если
# fixture-fingerprint gate прошёл (см. AMTICSEU_FINGERPRINT_LINES_KEV).
CALIBRATION_RANK_START: int = 15

# F-453-FU v2: peak'и-маркеры калибрационной фикстуры (AmTiCsEu/AmCs/AmCsEu).
# Calibration-tier активен IFF в peak list одновременно видны peak'и возле
# ВСЕХ линий ниже (Cs-137 661.66 + Am-241 59.54). На Ra-226/Th-232/K-40/Co-60
# нет ни Cs ни Am → gate fail → calibration-tier silent.
AMTICSEU_FINGERPRINT_LINES_KEV: Tuple[float, ...] = (661.66, 59.54)


# ──────────────────────────────────────────────────────────────────
# Express patterns (multi-line signatures)
# ──────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class ExpressPattern:
    name: str                                # short id
    nuclide: str                             # what it confirms
    required_lines_keV: Tuple[float, ...]    # all must be matched
    minimum_required: int                    # at least this many
    description: str


EXPRESS_PATTERNS: Tuple[ExpressPattern, ...] = (
    ExpressPattern(
        name="Co-60 pair",
        nuclide="Co-60",
        required_lines_keV=(1173.23, 1332.49),
        minimum_required=2,
        description="Co-60 doublet — оба пика должны быть видны вместе"
    ),
    ExpressPattern(
        name="Bi-214 Ra-chain pair",
        nuclide="Bi-214",
        required_lines_keV=(609.31, 1764.49),
        minimum_required=2,
        description="Bi-214 / радоновый ряд — пара 609 + 1764"
    ),
    ExpressPattern(
        name="Bi-214 quartet",
        nuclide="Bi-214",
        required_lines_keV=(609.31, 1120.29, 1764.49, 2204.21),
        minimum_required=3,
        description="Bi-214 quartet — ≥3 из 4 линий"
    ),
    ExpressPattern(
        name="Pb-214 doublet",
        nuclide="Pb-214",
        required_lines_keV=(351.93, 295.22),
        minimum_required=2,
        description="Pb-214 — обе главные линии 352 + 295"
    ),
    ExpressPattern(
        name="Th-232 strong",
        nuclide="Th-232 chain",
        required_lines_keV=(2614.51, 911.20),
        minimum_required=2,
        description="Минимальная подпись Th-232 — Tl-208 2615 + Ac-228 911"
    ),
    ExpressPattern(
        name="Th-232 triplet",
        nuclide="Th-232 chain",
        required_lines_keV=(2614.51, 911.20, 583.19),
        minimum_required=3,
        description="Полная подпись Th-232 — все 3 главные линии"
    ),
    ExpressPattern(
        name="Cs-137 + Ba IC",
        nuclide="Cs-137",
        required_lines_keV=(661.66, 32.07),
        minimum_required=2,
        description="Cs-137 + Ba Kα от внутр. конверсии"
    ),
)


# ──────────────────────────────────────────────────────────────────
# Result dataclasses
# ──────────────────────────────────────────────────────────────────

@dataclass
class AnchorMatch:
    """One found peak matched to one anchor entry."""
    anchor: AnchorEntry
    peak_channel: int
    peak_E_keV: float
    delta_keV: float
    sigma: float
    partner_required_but_missing: bool = False


@dataclass
class PatternConfirmation:
    pattern: ExpressPattern
    matched_lines_keV: List[float]
    missing_lines_keV: List[float]
    confirmed: bool                          # >= minimum_required
    note: str = ""


# ──────────────────────────────────────────────────────────────────
# Public matchers
# ──────────────────────────────────────────────────────────────────

def _amticseu_fingerprint_present(
    peak_E,
    fwhm_provider_keV,
    window_fwhm_multiple: float,
) -> bool:
    """F-453-FU v2 fixture-fingerprint gate.

    True ⇔ в peak list одновременно видны peak'и возле ВСЕХ линий из
    `AMTICSEU_FINGERPRINT_LINES_KEV` (Cs-137 661.66 + Am-241 59.54) —
    характерная подпись калибрационных фикстур AmTiCsEu/AmCs/AmCsEu.
    На Ra-226/Th-232/K-40/Co-60/природный фон — gate fail (нет ни Cs ни Am).

    Args:
        peak_E: список (peak, E_keV) от уже-калиброванного спектра
        fwhm_provider_keV: callable(E) -> FWHM в keV
        window_fwhm_multiple: half-width окна поиска (тот же, что для anchor)
    """
    needed = list(AMTICSEU_FINGERPRINT_LINES_KEV)
    found = [False] * len(needed)
    for _, e in peak_E:
        for i, line_E in enumerate(needed):
            if found[i]:
                continue
            tol = fwhm_provider_keV(line_E) * window_fwhm_multiple
            if abs(e - line_E) <= tol:
                found[i] = True
        if all(found):
            return True
    return all(found)


def find_anchor_matches(
    found_peaks,
    spec,
    *,
    fwhm_provider_keV,
    max_rank: int = 14,
    window_fwhm_multiple: float = 0.5,
) -> List[AnchorMatch]:
    """
    Try to match found peaks against the practical-anchor table.

    Returns one AnchorMatch per (anchor, peak) hit. A single peak can
    match multiple anchors (e.g. 240 keV peak hits both Pb-212 239 and
    Pb-214 242); in that case both matches are returned and the caller
    resolves the conflict via disambiguate.

    Args:
        found_peaks: iterable of FoundPeak from Mariscotti
        spec: Spectrum
        fwhm_provider_keV: callable(E) -> FWHM in keV
        max_rank: only consider anchors with rank <= this (default 14 = all)
        window_fwhm_multiple: match window half-width as multiple of FWHM
    """
    matches: List[AnchorMatch] = []
    # Index peaks by energy for fast lookup
    peak_E = [(p, spec.channel_to_energy(p.channel)) for p in found_peaks]
    # F-453-FU v2: fixture-fingerprint gate. Calibration-tier (rank >=
    # CALIBRATION_RANK_START) пробуется ТОЛЬКО на фикстурах с подписью
    # AmTiCsEu/AmCs — peak'и возле Cs-137 661.66 keV И Am-241 59.54 keV.
    # На Ra-226/Th-232/K-40/Co-60 — gate fail, calibration-tier silent.
    cal_tier_allowed = _amticseu_fingerprint_present(
        peak_E, fwhm_provider_keV, window_fwhm_multiple
    )
    for anchor in ANCHOR_RANKS:
        if anchor.rank > max_rank and anchor.rank < CALIBRATION_RANK_START:
            continue
        if anchor.rank >= CALIBRATION_RANK_START and not cal_tier_allowed:
            continue
        fwhm = max(2.0, fwhm_provider_keV(anchor.energy_keV))
        tol = fwhm * window_fwhm_multiple
        best: Optional[Tuple] = None
        for p, e in peak_E:
            d = abs(e - anchor.energy_keV)
            if d <= tol and (best is None or d < best[2]):
                best = (p, e, d)
        if best is not None:
            p, e, d = best
            # Partner check (Co-60-like pair anchors)
            partner_missing = False
            if anchor.requires_partner:
                partner_missing = True
                for part_E in anchor.partner_energies:
                    f2 = fwhm_provider_keV(part_E)
                    tol2 = f2 * window_fwhm_multiple
                    for _, e2 in peak_E:
                        if abs(e2 - part_E) <= tol2:
                            partner_missing = False
                            break
                    if not partner_missing:
                        break
            matches.append(AnchorMatch(
                anchor=anchor,
                peak_channel=p.channel,
                peak_E_keV=e,
                delta_keV=d,
                sigma=float(p.significance),
                partner_required_but_missing=partner_missing,
            ))
    return matches


def confirm_express_patterns(
    found_peaks,
    spec,
    *,
    fwhm_provider_keV,
    window_fwhm_multiple: float = 0.5,
) -> List[PatternConfirmation]:
    """
    For each EXPRESS_PATTERN, check how many of its required lines are
    represented in the found peaks.

    Returns one PatternConfirmation per pattern (whether confirmed or not).
    """
    peak_E = [(p, spec.channel_to_energy(p.channel)) for p in found_peaks]
    out: List[PatternConfirmation] = []
    for pat in EXPRESS_PATTERNS:
        matched: List[float] = []
        missing: List[float] = []
        for line_E in pat.required_lines_keV:
            tol = fwhm_provider_keV(line_E) * window_fwhm_multiple
            found = False
            for _, e in peak_E:
                if abs(e - line_E) <= tol:
                    found = True
                    break
            (matched if found else missing).append(line_E)
        confirmed = len(matched) >= pat.minimum_required
        note = ""
        if confirmed:
            note = (f"подтверждён: {len(matched)}/{len(pat.required_lines_keV)} линий "
                    f"({', '.join(f'{e:.0f}' for e in matched)})")
        else:
            note = (f"НЕ подтверждён: только {len(matched)}/{len(pat.required_lines_keV)} линий "
                    f"(missing: {', '.join(f'{e:.0f}' for e in missing)})")
        out.append(PatternConfirmation(
            pattern=pat,
            matched_lines_keV=matched,
            missing_lines_keV=missing,
            confirmed=confirmed,
            note=note,
        ))
    return out


def best_anchor_for_nuclide(nuclide: str) -> Optional[AnchorEntry]:
    """Highest-rank anchor for `nuclide`, or None."""
    cands = [a for a in ANCHOR_RANKS if a.nuclide == nuclide]
    if not cands:
        return None
    return min(cands, key=lambda a: a.rank)


def anchors_for_chain(chain_name: str) -> List[AnchorEntry]:
    """All anchors belonging to a given chain (e.g. 'Th-232' or 'U-238')."""
    return [a for a in ANCHOR_RANKS if a.chain == chain_name]


# ──────────────────────────────────────────────────────────────────
# F-87 — Anchor seeding bundle (Step 5α of SKILL.md workflow)
# ──────────────────────────────────────────────────────────────────
#
# Per the methodological refactor in v1.15.0, F-79 (anchor-rank) +
# F-80 (express-patterns) are properly part of Step 5 — they seed the
# calibration-verification (and optionally the calibration refit) with
# practical-anchor identifications. They were historically run inside
# Stage-1 identification; the call site in `staged_pipeline.py` now
# routes through this single helper so the position in the workflow
# matches the SKILL.md ordering.
#
# F-81 (7-line ЕРН check, Lsrm §9) remains a separate Pass C — it is
# the canonical *final calibration verification* on background- or ЕРН-
# rich spectra, not a seeding step.

@dataclass
class AnchorSeedResult:
    """Bundle of F-79 + F-80 outputs returned by ``seed_calibration_anchors``.

    Attributes
    ----------
    mode : str
        ``"sample"`` or ``"background"`` — passed through unchanged so
        callers can branch downstream without re-deriving the mode.
    anchor_matches : list of AnchorMatch
        F-79 hits against the practical-visibility table.
    pattern_confirmations : list of PatternConfirmation
        F-80 express-pattern verdicts (one per pattern; ``confirmed``
        flag tells whether the pattern reached its ``minimum_required``).
    """
    mode: str
    anchor_matches: List[AnchorMatch] = field(default_factory=list)
    pattern_confirmations: List[PatternConfirmation] = field(
        default_factory=list
    )


def seed_calibration_anchors(
    found_peaks,
    spec,
    *,
    mode: str = "sample",
    fwhm_provider_keV,
    max_rank: Optional[int] = None,
    window_fwhm_multiple: float = 0.5,
) -> AnchorSeedResult:
    """Run the F-79 + F-80 anchor-seeding bundle (Step 5α).

    This is the single entry point for the express identification used
    to seed Step 5 calibration verification / refit. The function is a
    thin orchestrator on top of :func:`find_anchor_matches` and
    :func:`confirm_express_patterns` — there is no new algorithmic
    work, only a consistent boundary so the call site in
    ``analyze_lsrm_spe`` does not duplicate logic and so test code can
    invoke the heuristic standalone.

    Parameters
    ----------
    found_peaks : iterable of FoundPeak
        Peaks from the channel-space Mariscotti search (Step 3).
    spec : Spectrum
        Provides ``channel_to_energy`` for window matching. The current
        stored / bootstrap energy calibration is used as-is; the
        recalibration hook (F-87c/d) is the caller's responsibility.
    mode : str, default ``"sample"``
        ``"sample"`` runs both F-79 (anchors) and F-80 (patterns) at
        ``max_rank = 10`` (skips 511 / U-235—Ra-226 ambiguity).
        ``"background"`` still runs F-79 + F-80 for ЕРН anchors but
        extends ``max_rank`` to 12 because background spectra rely on
        Pb-212 238.6 and Tl-208 583.2 as secondary anchors. The F-81
        7-line ЕРН check stays separate as the final calibration
        verification (Lsrm §9) — that call is **not** subsumed here.
    fwhm_provider_keV : callable
        ``f(E_keV) -> FWHM_keV``.
    max_rank : int, optional
        Override the mode-derived rank cap.
    window_fwhm_multiple : float, default 0.5
        Match window half-width as multiple of FWHM (passed through to
        both inner functions).

    Returns
    -------
    AnchorSeedResult
    """
    if mode not in {"sample", "background"}:
        raise ValueError(
            f"seed_calibration_anchors: mode must be "
            f"'sample' or 'background', got {mode!r}"
        )

    if max_rank is None:
        max_rank = 10 if mode == "sample" else 12

    anchor_matches = find_anchor_matches(
        found_peaks, spec,
        fwhm_provider_keV=fwhm_provider_keV,
        max_rank=max_rank,
        window_fwhm_multiple=window_fwhm_multiple,
    )

    pattern_confirmations = confirm_express_patterns(
        found_peaks, spec,
        fwhm_provider_keV=fwhm_provider_keV,
        window_fwhm_multiple=window_fwhm_multiple,
    )

    return AnchorSeedResult(
        mode=mode,
        anchor_matches=anchor_matches,
        pattern_confirmations=pattern_confirmations,
    )


# ──────────────────────────────────────────────────────────────────
# F-88 — User-prioritized express order + chain dominance
# ──────────────────────────────────────────────────────────────────
#
# Per user methodology (2026-05-29):
#
#   "Важен порядок эвристики при экспресс идентификации:
#    1. 2615 → торий. Это фиксирует жёстко его наличие; уже этого
#       достаточно для первичной калибровки. Как правило остальные
#       изотопы при выраженном тории идентифицировать на этапе
#       калибровки бессмысленно. Данные о наличии тория должны
#       жёстко передаваться на этап идентификации пиков.
#    2. 1461 → калий. При наличии тория возможно наложение [с Ac-228 1459].
#    3.  662 → цезий
#    4. 1173+1332 → кобальт
#    5.  609+1764 → радон/урановый ряд
#    6.   59.5  → америций"
#
# This priority differs subtly from ANCHOR_RANKS (which is sorted by
# *practical visibility on NaI*) — the user's order is by **express
# diagnostic value at calibration time**. Th is the trump card; once
# 2615 is found, the chain becomes a hard prior for downstream
# identification. K-40 must then be confirmed against the Ac-228 1459
# overlap risk.

@dataclass(frozen=True)
class PrioritySignal:
    """One entry in the user-ordered express anchor priority list."""
    order: int                       # 1..6 (1 = highest priority)
    label: str                       # short Russian label for reports
    nuclide_or_chain: str            # what this signal identifies
    chain: str                       # "Th-232" / "U-238" / ""
    required_lines_keV: Tuple[float, ...]
    minimum_required: int            # how many required_lines must hit
    note: str                        # rationale in English


# The canonical user-supplied order. Used by derive_priority_findings.
USER_PRIORITY_ORDER: Tuple[PrioritySignal, ...] = (
    PrioritySignal(
        order=1,
        label="Th-232 (Tl-208 2615 keV)",
        nuclide_or_chain="Tl-208",
        chain="Th-232",
        required_lines_keV=(2614.51,),
        minimum_required=1,
        note="Trump card — locks Th-232 chain dominance for the rest "
             "of the workflow. Hard-passes to Step 7 identification.",
    ),
    PrioritySignal(
        order=2,
        label="K-40 1461 keV",
        nuclide_or_chain="K-40",
        chain="",
        required_lines_keV=(1460.82,),
        minimum_required=1,
        note="Single ⁴⁰K line. When Th-232 dominant, overlap with "
             "Ac-228 1459.20 keV (I=0.85%) must be flagged because "
             "NaI 63×63 cannot resolve the doublet.",
    ),
    PrioritySignal(
        order=3,
        label="Cs-137 662 keV",
        nuclide_or_chain="Cs-137",
        chain="",
        required_lines_keV=(661.66,),
        minimum_required=1,
        note="Bright single-line FEP in a clean region.",
    ),
    PrioritySignal(
        order=4,
        label="Co-60 1173+1332 keV (pair)",
        nuclide_or_chain="Co-60",
        chain="",
        required_lines_keV=(1173.23, 1332.49),
        minimum_required=2,
        note="Paired signature — both lines must be present.",
    ),
    PrioritySignal(
        order=5,
        label="Ra/U chain (Bi-214 609+1764 keV)",
        nuclide_or_chain="Bi-214",
        chain="U-238",
        required_lines_keV=(609.31, 1764.49),
        minimum_required=2,
        note="Bi-214 Ra-chain pair confirms the U-238 / Rn-222 ingrowth.",
    ),
    PrioritySignal(
        order=6,
        label="Am-241 59.5 keV",
        nuclide_or_chain="Am-241",
        chain="",
        required_lines_keV=(59.54,),
        minimum_required=1,
        note="Characteristic low-energy line.",
    ),
)


@dataclass
class PriorityFinding:
    """Per-signal outcome of derive_priority_findings."""
    signal: PrioritySignal
    matched: bool                              # >= minimum_required hits
    matched_lines_keV: Tuple[float, ...] = ()
    missing_lines_keV: Tuple[float, ...] = ()
    significance: float = 0.0                  # max σ across matched lines
    note: str = ""


def derive_priority_findings(
    anchor_matches: List[AnchorMatch],
    pattern_confirmations: List[PatternConfirmation],
) -> List[PriorityFinding]:
    """Evaluate USER_PRIORITY_ORDER against the F-79 + F-80 outputs.

    Returns one PriorityFinding per entry in USER_PRIORITY_ORDER, in
    the same order (so the report can list them top-to-bottom by
    diagnostic value).

    The matching logic:
      * For each required line, check if any AnchorMatch has
        ``peak_E_keV`` within 1 keV of the library energy (the
        AnchorMatch already passed its FWHM-based window test).
      * If the signal is part of an EXPRESS_PATTERN (Co-60 pair,
        Bi-214 pair), use the pattern confirmation as the authority —
        this matches the F-80 semantics exactly.
    """
    # Build a lookup of (E_lib, peak_E, σ) from anchor matches
    anchor_hits: List[Tuple[float, float, float]] = []
    for am in anchor_matches:
        if not am.anchor.nuclide:
            # Skip annihilation / ambiguous anchors
            continue
        anchor_hits.append(
            (am.anchor.energy_keV, am.peak_E_keV, am.sigma)
        )

    # Build a lookup of pattern confirmations by nuclide
    pat_by_nuclide = {}
    for pc in pattern_confirmations:
        # Use the most permissive confirmation per nuclide
        existing = pat_by_nuclide.get(pc.pattern.nuclide)
        if existing is None or (pc.confirmed and not existing.confirmed):
            pat_by_nuclide[pc.pattern.nuclide] = pc

    findings: List[PriorityFinding] = []
    for sig in USER_PRIORITY_ORDER:
        # Try express-pattern path first (Co-60, Bi-214 pair, Th-strong)
        pattern_path_match = None
        if sig.minimum_required >= 2:
            # Try matching pattern_confirmations for this nuclide/chain
            pat = pat_by_nuclide.get(sig.nuclide_or_chain)
            if pat is not None and pat.confirmed:
                pattern_path_match = pat

        matched_lines: List[float] = []
        missing_lines: List[float] = []
        max_sigma = 0.0
        for line_E in sig.required_lines_keV:
            hit = False
            for E_lib, peak_E, sigma in anchor_hits:
                # Allow 1 keV tolerance on E_lib match (anchor was
                # already FWHM-checked) AND the peak must be close.
                if abs(E_lib - line_E) <= 1.0:
                    hit = True
                    if sigma > max_sigma:
                        max_sigma = sigma
                    break
            if hit:
                matched_lines.append(line_E)
            else:
                missing_lines.append(line_E)

        matched = (
            len(matched_lines) >= sig.minimum_required
            or pattern_path_match is not None
        )

        if matched:
            note = (
                f"MATCHED ({len(matched_lines)}/"
                f"{len(sig.required_lines_keV)} lines, "
                f"σ≤{max_sigma:.1f})"
            )
        else:
            note = (
                f"missing: {', '.join(f'{e:.0f}' for e in missing_lines)}"
            )

        findings.append(PriorityFinding(
            signal=sig,
            matched=matched,
            matched_lines_keV=tuple(matched_lines),
            missing_lines_keV=tuple(missing_lines),
            significance=max_sigma,
            note=note,
        ))

    return findings


@dataclass
class ChainDominance:
    """Outcome of derive_chain_dominance — what nuclide chain controls
    the spectrum, and the evidence used to decide.

    A ``True`` flag is a **hard prior** for Step 7 identification:
    chain members must be considered confirmed candidates even when
    individual χ² tests are tight.

    F-89d / v1.15.2 — when the filename binds the source to a specific
    chain (e.g. ``Th232_*.spe``), competing chains require stricter
    evidence to be marked dominant. See ``suppressed_chains`` and
    ``suppression_reason`` for transparency.
    """
    th232: bool = False
    u238: bool = False
    th232_evidence: Tuple[str, ...] = ()
    u238_evidence: Tuple[str, ...] = ()
    th232_strength_sigma: float = 0.0
    u238_strength_sigma: float = 0.0
    reason: str = ""
    # F-89d — chains that were *evidence-positive* but suppressed by
    # the filename-binding rule. Used by the orchestrator to drop
    # spurious proxies from final_detected.
    suppressed_chains: Tuple[str, ...] = ()
    suppression_reason: str = ""


# Th-232 chain members to test (anchor energies present in ANCHOR_RANKS)
_TH232_CHAIN_LINES = {
    2614.51: "Tl-208 2614.51 (rank 1)",
    911.20:  "Ac-228 911.20 (rank 7)",
    583.19:  "Tl-208 583.19 (rank 11)",
    238.63:  "Pb-212 238.63 (rank 12)",
}

_U238_CHAIN_LINES = {
    1764.49: "Bi-214 1764.49 (rank 6)",
    609.31:  "Bi-214 609.31 (rank 8)",
    351.93:  "Pb-214 351.93 (rank 9)",
    295.22:  "Pb-214 295.22 (rank 10)",
}


def derive_chain_dominance(
    anchor_matches: List[AnchorMatch],
    pattern_confirmations: List[PatternConfirmation],
    *,
    filename_chains_claimed: Optional[set] = None,
) -> ChainDominance:
    """Decide whether Th-232 and/or U-238 chains dominate the spectrum.

    Rules:

    Th-232 dominant when **any** of:
      * Tl-208 2614.51 matched with σ ≥ 5  (trump card per user
        methodology — 2615 alone is enough to fix Th presence)
      * ≥ 2 distinct Th-chain anchors matched (Tl-208 2615, Ac-228
        911, Tl-208 583, Pb-212 238)
      * "Th-232 strong" or "Th-232 triplet" express pattern confirmed

    U-238 dominant when **any** of:
      * "Bi-214 Ra-chain pair" express pattern confirmed (609+1764)
      * "Bi-214 quartet" express pattern confirmed (≥3 of 609, 1120,
        1764, 2204)
      * ≥ 3 distinct U-chain anchors matched

    Both can be True simultaneously (natural background often has both
    Th and U chains active).

    F-89d / v1.15.2 — filename-binding suppression rule.
    When ``filename_chains_claimed`` indicates the source is bound to
    a specific chain set (e.g. ``{"Th-232"}`` from a ``Th232_*.spe``
    filename) AND a competing chain is *not* explicitly claimed, the
    stricter evidence threshold applies:

      * On a Th-only filename, U-238 dominance requires the Bi-214
        **quartet** (≥3 of 609, 1120, 1764, 2204 matched). The
        Ra-pair (609+1764 alone) is insufficient because on NaI 63×63
        the 609 keV peak is often Tl-208 583 keV shifted by Compton
        overlap, and 1764 may be a Bi-212 1620 keV tail.
      * Same symmetric rule for U-only filenames vs Th: Th-232 requires
        ≥3 distinct Th-chain anchors (trump-card alone insufficient).

    Suppressed chain identity is recorded in ``suppressed_chains`` and
    ``suppression_reason`` so the orchestrator can drop the
    corresponding chain proxies from ``final_detected``.
    """
    # Count Th-232 evidence
    th_evidence: List[str] = []
    th_sigma = 0.0
    th_2614_sigma = 0.0
    th_matched_lines = set()
    for am in anchor_matches:
        if am.anchor.chain == "Th-232":
            E = am.anchor.energy_keV
            if E in _TH232_CHAIN_LINES and E not in th_matched_lines:
                th_matched_lines.add(E)
                th_evidence.append(
                    f"{_TH232_CHAIN_LINES[E]}, σ={am.sigma:.1f}"
                )
                if am.sigma > th_sigma:
                    th_sigma = am.sigma
                if abs(E - 2614.51) < 1.0 and am.sigma > th_2614_sigma:
                    th_2614_sigma = am.sigma

    # Th-232 dominance decision
    th_dom = False
    th_reason = ""
    if th_2614_sigma >= 5.0:
        th_dom = True
        th_reason = (
            f"trump-card rule: Tl-208 2614.51 matched at σ={th_2614_sigma:.1f}"
        )
    elif len(th_matched_lines) >= 2:
        th_dom = True
        th_reason = (
            f"multi-anchor rule: {len(th_matched_lines)} Th-chain "
            f"anchors matched"
        )
    else:
        # Try express-pattern path
        for pc in pattern_confirmations:
            if pc.confirmed and pc.pattern.nuclide == "Th-232 chain":
                th_dom = True
                th_reason = (
                    f"express-pattern rule: '{pc.pattern.name}' confirmed"
                )
                break

    # Count U-238 evidence
    u_evidence: List[str] = []
    u_sigma = 0.0
    u_matched_lines = set()
    for am in anchor_matches:
        if am.anchor.chain == "U-238":
            E = am.anchor.energy_keV
            if E in _U238_CHAIN_LINES and E not in u_matched_lines:
                u_matched_lines.add(E)
                u_evidence.append(
                    f"{_U238_CHAIN_LINES[E]}, σ={am.sigma:.1f}"
                )
                if am.sigma > u_sigma:
                    u_sigma = am.sigma

    u_dom = False
    u_reason = ""
    # Check Bi-214 express patterns first (user-prioritized signal #5)
    bi_pair_confirmed = False
    bi_quartet_confirmed = False
    for pc in pattern_confirmations:
        if not pc.confirmed:
            continue
        if pc.pattern.name == "Bi-214 Ra-chain pair":
            bi_pair_confirmed = True
        elif pc.pattern.name == "Bi-214 quartet":
            bi_quartet_confirmed = True

    if bi_pair_confirmed or bi_quartet_confirmed:
        u_dom = True
        u_reason = (
            f"express-pattern rule: '"
            f"{'Bi-214 quartet' if bi_quartet_confirmed else 'Bi-214 Ra-chain pair'}"
            f"' confirmed"
        )
    elif len(u_matched_lines) >= 3:
        u_dom = True
        u_reason = (
            f"multi-anchor rule: {len(u_matched_lines)} U-chain "
            f"anchors matched"
        )

    # ─── F-89d filename-binding suppression ──────────────────────
    # Per user methodology (2026-05-29): "Радия в образце точно нет.
    # Его следы могут быть только от фона, но cps фона пренебрежимо
    # мал." On a Th-only source, U-chain identification is presumed
    # to be NaI Compton-confusion (Tl-208 583 ≈ Bi-214 609 within
    # one FWHM on Gamma-1S), NOT real U/Ra contamination. The rule
    # fires UNCONDITIONALLY when the filename excludes the chain —
    # not just when the chain dominance would have fired.
    suppressed: List[str] = []
    suppression_reason = ""
    if filename_chains_claimed:
        th_in_filename = "Th-232" in filename_chains_claimed
        u_in_filename = (
            "U-238" in filename_chains_claimed
            or "Ra-226" in filename_chains_claimed
        )
        # Single-isotope sources (Cs-137, K-40, Co-60, Am-241) bind
        # the source to that isotope only — neither Th nor U chains
        # are claimed. The chain-suppression rules don't fire on
        # those (the filename hint just drives candidate selection).
        # Suppression activates only when the filename CLAIMS one
        # chain explicitly.
        # Th-only filename: U-238 chain suppressed unless quartet.
        if th_in_filename and not u_in_filename and not bi_quartet_confirmed:
            suppressed.append("U-238")
            suppression_reason = (
                "Filename binds source to Th-232 chain only. U-238 evidence "
                "(Bi-214 Ra-pair 609+1764) is suppressed because on NaI 63×63 "
                "the 609 keV peak is often Tl-208 583 keV shifted by Compton "
                "overlap. U-238 dominance would require the Bi-214 quartet "
                "(≥3 of 609, 1120, 1764, 2204) — not satisfied here."
            )
            u_dom = False
            u_evidence = []
            u_reason = "suppressed by filename-binding rule"
        # U-only filename: Th-232 requires multi-anchor (not just trump)
        if u_in_filename and not th_in_filename and len(th_matched_lines) < 3:
            suppressed.append("Th-232")
            suppression_reason = (
                "Filename binds source to U-238 chain only. Th-232 evidence "
                "(single trump-card anchor) is suppressed; full Th-232 "
                "dominance would require ≥3 distinct Th-chain anchors."
            )
            th_dom = False
            th_evidence = []
            th_reason = "suppressed by filename-binding rule"

    # Combined reason for the dominance verdict
    full_reason = []
    if th_dom:
        full_reason.append(f"Th-232: {th_reason}")
    if u_dom:
        full_reason.append(f"U-238: {u_reason}")
    if not th_dom and not u_dom and not suppressed:
        full_reason.append("neither chain dominant")
    if suppressed:
        full_reason.append(f"suppressed: {'+'.join(suppressed)}")

    return ChainDominance(
        th232=th_dom,
        u238=u_dom,
        th232_evidence=tuple(th_evidence),
        u238_evidence=tuple(u_evidence),
        th232_strength_sigma=th_sigma,
        u238_strength_sigma=u_sigma,
        reason=" | ".join(full_reason),
        suppressed_chains=tuple(suppressed),
        suppression_reason=suppression_reason,
    )


# Chain proxy nuclides — when chain dominance is set, these become
# strong-prior candidates at Step 7 even if normal ЕРН filtering
# wouldn't include them.
TH232_PROXY_NUCLIDES: Tuple[str, ...] = (
    "Tl-208", "Pb-212", "Ac-228", "Bi-212", "Pb-214",
)
U238_PROXY_NUCLIDES: Tuple[str, ...] = (
    "Bi-214", "Pb-214", "Pb-210", "Ra-226",
)


__all__ = [
    "AnchorEntry", "ANCHOR_RANKS",
    "ExpressPattern", "EXPRESS_PATTERNS",
    "AnchorMatch", "PatternConfirmation",
    "find_anchor_matches", "confirm_express_patterns",
    "best_anchor_for_nuclide", "anchors_for_chain",
    "AnchorSeedResult", "seed_calibration_anchors",
    # F-88
    "PrioritySignal", "USER_PRIORITY_ORDER",
    "PriorityFinding", "derive_priority_findings",
    "ChainDominance", "derive_chain_dominance",
    "TH232_PROXY_NUCLIDES", "U238_PROXY_NUCLIDES",
]
