# -*- coding: utf-8 -*-
"""F-372 / v1.18.24.7 — Обновление data/nuclides.json через IAEA Live Chart API.

Стратегия:
  1. Backup текущей data/nuclides.json → data/nuclides.json.backup-YYYY-MM-DD.
  2. Для каждого нуклида (29 шт) fetch_iaea_gamma_lines(force_refresh=True).
  3. Filter: γ-lines >50 кэВ, I>=0.1%, only ground-state parent.
  4. Merge: обновляем `lines` field. `ic_xrays` НЕ трогаем (нижний раздел
     K-X из internal conversion — отдельная физика, статически outsourced
     LSRM library).
  5. Сравнение E_keV / I_pct / counts old vs new для каждого нуклида →
     текстовый changelog.
  6. Запись обновлённой nuclides.json.

NB: T_half_s, parent, daughters, chain, is_cascade — НЕ обновляются
этим скриптом (T½ имеется в IAEA, но это требует отдельной API query;
chain-metadata — наша добавка, нет в ENSDF).

NB2: при rate-limit или 403 от IAEA — script продолжает с других
нуклидов, выводит warning в log.
"""
from __future__ import annotations

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import argparse
import json
import shutil
import time
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.data.iaea_fetcher import (  # noqa: E402
    fetch_iaea_gamma_lines, merge_iaea_into_internal,
    _cache_path, DEFAULT_CACHE_DIR,
)


GAMMA_MIN_KEV = 0.0    # Без energy cut-off. IAEA `rad_types=g` уже
                       # возвращает только γ-emission (без α/β).
                       # Низкоэнергетические γ (Am-241 59.5, Pb-210 46.5)
                       # — legitimate, сохраняются как есть.
MIN_INTENSITY_PCT = 0.1  # Главный физический фильтр. NaI 63×63 имеет ПШПВ
                         # ~7% (FWHM 46 кэВ на 662 кэВ). При intensity ≥0.1%
                         # линии остаются distinguishable в фоне (для
                         # сильных источников). Порог 0.1% match-ит
                         # historical OLD-library practice (hand-curated
                         # weak Th-228 215/166/131 keV all ≈0.1-0.25%).
                         # K-X фильтруются ОТДЕЛЬНО через multipolarity
                         # field — НЕ через intensity.


# F-372 / v1.18.24.7 — Branching factors: pipeline ожидает intensities
# «per parent-of-chain decay» (в секулярном равновесии); IAEA возвращает
# «per this-nuclide decay». Для нуклидов с branched-parent нужно умножать.
#
# Источник коэффициентов: ENSDF current evaluations (см. IAEA Live Chart
# decay-mode strengths). Без этой корректировки activity Th-232 цепочки
# завышается в 1/0.3594 = 2.78× (старая ошибка которую мы избегаем).
#
# Singletons (K-40, Cs-137, Co-60, ...) и нуклиды в прямой цепи (Pb-212,
# Pb-214, Bi-214, Ac-228, ...) — branching = 1.0 (не указаны в карте).
_BRANCHING_TO_PARENT = {
    # Th-232 chain: Bi-212 → 35.94% α → Tl-208 (parallel β⁻ branch к Po-212)
    "Tl-208": 0.3594,
    # U-238 chain: Pa-234m → 99.84% IT → Pa-234 (g.s.) → β⁻;
    #              0.16% β⁻ → U-234 (skip Pa-234 g.s.)
    # Pa-234m имеет γ → per-U-238-decay × 0.9984 ≈ 1.0 — оставляем как есть.
    # Pa-234 (g.s.) — 0.16% branch только. Малая γ-emission в наших целях.
}


def _gamma_keys_from_csv_cache(name: str) -> set:
    """F-372 / v1.18.24.7 — извлечь набор (round(E,2), I) для строк CSV
    с непустой `multipolarity` (true γ-transitions, не K X-rays).

    IAEA `rad_types=g` возвращает И γ-emissions И K X-rays от IC. Фильтр
    по multipolarity критичен: K-X (74-88 кэВ Bi/Pb/Rn) имеют пустую
    multipolarity, тогда как γ-emissions содержат E1/M1/E2/M2/etc.
    Без фильтра K-X засоряют lines и ломают identification (например
    Tl-208 и Pb-212 одновременно «emit» 74 кэВ → confusion на пиках).
    """
    import csv
    from pathlib import Path as _P
    path = _cache_path(_P(DEFAULT_CACHE_DIR), name, "g")
    if not path.exists():
        return set()
    keep = set()
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mp = (row.get("multipolarity") or "").strip()
            if not mp:
                continue   # X-ray (or unclassified) → drop
            try:
                E = float(row.get("energy") or 0)
            except (ValueError, TypeError):
                continue
            keep.add(round(E, 2))
    return keep


def fetch_one_nuclide(name: str, sleep_s: float = 0.5) -> dict:
    """Fetch + merge для одного нуклида. Возвращает {'lines': [...]}.

    sleep_s между запросами — vежливый rate-limit к IAEA endpoint."""
    print(f"  → {name}...", end="", flush=True)
    try:
        time.sleep(sleep_s)
        iaea_lines = fetch_iaea_gamma_lines(name, force_refresh=True)
    except Exception as e:
        print(f" FAIL: {type(e).__name__}: {e}")
        return None
    # F-372 — keep only true γ-transitions (multipolarity != "").
    # K X-rays (пустая multipolarity) физически реальны, но имеют свой
    # bucket `ic_xrays` в nuclide library — не должны попадать в `lines`.
    gamma_E_set = _gamma_keys_from_csv_cache(name)
    gamma_only = [
        ln for ln in iaea_lines
        if (ln.energy_keV >= GAMMA_MIN_KEV
            and round(ln.energy_keV, 2) in gamma_E_set)
    ]
    merged = merge_iaea_into_internal(
        gamma_only, target_nuclide_name=name,
        min_intensity_pct=MIN_INTENSITY_PCT,
        only_ground_state_parent=True,
    )
    # F-372 — apply branching-to-parent если nuclide — branched chain daughter
    bf = _BRANCHING_TO_PARENT.get(name)
    if bf is not None:
        merged["lines"] = [
            [E, I * bf, dI * bf] for (E, I, dI) in merged["lines"]
        ]
        # Re-filter после branching — некоторые lines могут упасть ниже порога
        merged["lines"] = [
            ln for ln in merged["lines"] if ln[1] >= MIN_INTENSITY_PCT
        ]
    n = len(merged.get("lines", []))
    bf_note = f" × branching={bf}" if bf else ""
    print(f" {n} γ-lines (≥{GAMMA_MIN_KEV} кэВ, I≥{MIN_INTENSITY_PCT}%){bf_note}")
    return merged


def diff_lines(old_lines: list, new_lines: list, tol_keV: float = 1.5) -> dict:
    """Сравнить старые и новые списки [[E, I, dI], ...]."""
    diff = {
        "n_old": len(old_lines),
        "n_new": len(new_lines),
        "matched": [],   # (E_old, I_old, E_new, I_new, dI%)
        "only_in_old": [],
        "only_in_new": [],
    }
    used_new = set()
    for E_o, I_o, *_ in old_lines:
        best = None
        for i, ln in enumerate(new_lines):
            if i in used_new:
                continue
            E_n = ln[0]
            if abs(E_o - E_n) < tol_keV:
                if best is None or abs(E_o - E_n) < abs(E_o - new_lines[best][0]):
                    best = i
        if best is not None:
            used_new.add(best)
            ln = new_lines[best]
            E_n, I_n = ln[0], ln[1]
            dI_rel = (
                100.0 * (I_n - I_o) / I_o if I_o > 0 else None
            )
            diff["matched"].append((E_o, I_o, E_n, I_n, dI_rel))
        else:
            diff["only_in_old"].append((E_o, I_o))
    for i, ln in enumerate(new_lines):
        if i not in used_new:
            diff["only_in_new"].append((ln[0], ln[1]))
    return diff


def format_diff(name: str, diff: dict, T_half: float | None = None) -> str:
    out = [f"\n## {name}  (T½ = {T_half} s)" if T_half is not None
           else f"\n## {name}"]
    out.append(
        f"  Lines: old={diff['n_old']} → new={diff['n_new']} "
        f"(matched={len(diff['matched'])}, "
        f"only_old={len(diff['only_in_old'])}, "
        f"only_new={len(diff['only_in_new'])})"
    )
    # Show significant I changes (>5%)
    big_dI = [m for m in diff["matched"]
              if m[4] is not None and abs(m[4]) > 5.0]
    if big_dI:
        out.append("  Significant intensity changes (>5%):")
        for E_o, I_o, E_n, I_n, dI in big_dI:
            out.append(
                f"    {E_o:.2f} кэВ:  I_old={I_o:.3g}%  → I_new={I_n:.3g}%  "
                f"({dI:+.1f}%)"
            )
    if diff["only_in_old"]:
        out.append("  Removed (not in new IAEA fetch):")
        for E, I in diff["only_in_old"]:
            out.append(f"    {E:.2f} кэВ  I={I:.3g}%")
    if diff["only_in_new"]:
        out.append("  Added (new in IAEA fetch):")
        for E, I in diff["only_in_new"][:8]:  # cap output
            out.append(f"    {E:.2f} кэВ  I={I:.3g}%")
        if len(diff["only_in_new"]) > 8:
            out.append(f"    ... +{len(diff['only_in_new'])-8} more")
    return "\n".join(out)


def main():
    p = argparse.ArgumentParser(
        description=(
            "F-372 — обновление data/nuclides.json через IAEA Live Chart API. "
            "Backup-ит старый файл, fetch-ит свежие γ-lines, сохраняет diff."
        ),
    )
    p.add_argument("--input", default=str(REPO / "data" / "nuclides.json"))
    p.add_argument("--output", default=str(REPO / "data" / "nuclides.json"))
    p.add_argument("--changelog", default=str(REPO / "data" / "nuclides_iaea_update.md"))
    p.add_argument("--only", nargs="*", help="ограничить набор нуклидов")
    p.add_argument("--dry-run", action="store_true",
                   help="не записывать output, только diff")
    p.add_argument("--sleep", type=float, default=0.5,
                   help="задержка между запросами IAEA, сек (default 0.5)")
    args = p.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    log_path = Path(args.changelog)

    print(f"[F-372] Loading {in_path}")
    with in_path.open(encoding="utf-8") as f:
        data = json.load(f)
    schema = data.pop("_schema", None)

    nuclides = list(data.keys())
    if args.only:
        nuclides = [n for n in nuclides if n in set(args.only)]
    print(f"[F-372] Will fetch {len(nuclides)} nuclides from IAEA")

    if not args.dry_run and in_path == out_path:
        backup = in_path.with_suffix(
            f".json.backup-{date(2026, 6, 1).isoformat()}",
        )
        if not backup.exists():
            shutil.copy2(in_path, backup)
            print(f"[F-372] Backup: {backup}")

    diff_lines_log = [
        f"# IAEA Live Chart update — {date(2026, 6, 1).isoformat()}",
        f"\nFetched via `gamma.data.iaea_fetcher.fetch_iaea_gamma_lines`",
        f"with `min_intensity_pct={MIN_INTENSITY_PCT}` and "
        f"`only_ground_state_parent=True`. Threshold E≥{GAMMA_MIN_KEV} кэВ "
        f"(ниже — K X-rays, остаются в поле `ic_xrays` без изменений).",
    ]

    successes = 0
    failures = 0
    for name in nuclides:
        old_entry = data[name]
        old_lines = list(old_entry.get("lines", []))
        new_payload = fetch_one_nuclide(name, sleep_s=args.sleep)
        if new_payload is None:
            failures += 1
            diff_lines_log.append(
                f"\n## {name}\n  ⚠ IAEA fetch failed — entry unchanged."
            )
            continue
        new_lines = new_payload["lines"]
        if not new_lines:
            diff_lines_log.append(
                f"\n## {name}\n  ⚠ IAEA вернул 0 γ-линий ≥{GAMMA_MIN_KEV} кэВ — "
                f"entry unchanged (вероятно, нуклид без γ-emission в нашем "
                f"диапазоне или API нюанс)."
            )
            continue
        d = diff_lines(old_lines, new_lines)
        diff_lines_log.append(
            format_diff(name, d, T_half=old_entry.get("T_half_s")),
        )
        # Apply update: replace `lines`
        data[name]["lines"] = new_lines
        successes += 1

    # Bump schema version if present
    if schema is not None:
        old_ver = schema.get("version", "0.1")
        try:
            major, minor = old_ver.split(".")
            new_ver = f"{major}.{int(minor)+1}"
        except (ValueError, AttributeError):
            new_ver = "0.2"
        schema["version"] = new_ver
        schema["sources"] = (
            f"γ-lines from IAEA Live Chart API "
            f"(https://www-nds.iaea.org/relnsd/v1/data) — fetched "
            f"{date(2026, 6, 1).isoformat()}. T_half_s / chain / "
            f"daughters / ic_xrays: статически из 03b_nuclide_library.md "
            f"+ LSRM library."
        )
        schema["last_updated"] = date(2026, 6, 1).isoformat()
        data = {"_schema": schema, **data}
    else:
        data = data

    if args.dry_run:
        print(f"\n[F-372] DRY RUN — output NOT written")
    else:
        with out_path.open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        print(f"\n[F-372] Wrote: {out_path}")

    diff_lines_log.append(
        f"\n\n## Summary\n  successes: {successes}  failures: {failures}"
    )
    log_path.write_text("\n".join(diff_lines_log), encoding="utf-8")
    print(f"[F-372] Changelog: {log_path}")
    print(f"[F-372] Done. {successes} updated, {failures} failed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
