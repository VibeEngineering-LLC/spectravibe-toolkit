# -*- coding: utf-8 -*-
"""F-448 — раннер МЭС/template-метода + блок сравнения с efficiency-методом.

Запускает количественную оценку активности по методу эталонных спектров
(МЭС, LSRM §12): строит матрицу чувствительности R из эталонов Поверки
2024, решает линейную систему по пробе → A_k(МЭС) с σ. Опционально
сравнивает с efficiency/peak-area результатами (из JSON прогона
run_plan_a.py) и выдаёт таблицу A_eff vs A_МЭС vs паспорт (Δ%/|z|).

Использование (из корня репо):

    # Только МЭС по умолчательной пробе Th-232 420-7-17:
    PYTHONIOENCODING=utf-8 python scripts/run_template_assay.py

    # С compare против efficiency-прогона:
    PYTHONIOENCODING=utf-8 python scripts/run_template_assay.py \
        --efficiency-json demo_reports/<ts>_<label>/cert_zcheck.json \
        --out-dir demo_reports/<ts>_<label>

Env / флаги:
    --sample PATH         путь к пробе .spe (default: Th232_420-7-17)
    --bg PATH             путь к фону (default: Фон закр кр вода_13)
    --mass KG             масса пробы, кг (для записи в отчёт; default 1.6)
    --efficiency-json P   cert_zcheck.json ИЛИ <stem>_report.json прогона
                          run_plan_a.py — источник A_efficiency для compare
    --cert-A / --cert-rel паспорт пробы (default 1940 ±6% = Th-232 ОИСН-16)
    --out-dir DIR         куда писать template_compare.json/.md
                          (default: рядом с efficiency-json, иначе
                          demo_reports/<ts>_template_assay/)

ОГРАНИЧЕНИЯ (HARD): этот раннер НЕ генерит operator-facing HTML-отчёт.
Главный отчёт — только через run_plan_a.py. Здесь — отдельный
compare-артефакт (JSON + MD) для того же demo_reports-каталога.

Циркулярность: ассерт sample_path != Th-эталон (420-17031). Th-эталон —
независимый источник 860 Бк/кг, проба — 1940 Бк/кг.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))  # make `gamma` importable when not installed

from gamma.io.lsrm_spe import read_lsrm_spe
from gamma.activity.template_method import template_assay, BackgroundSpec
from gamma.activity.build_sensitivity_matrix import (
    MES_WINDOWS,
    windows_keV_only,
    SAMPLE_TH232_FILE,
    DEFAULT_BG_FILE,
    ETALON_FILES,
    matrix_descriptor,
    build_sensitivity_matrix_from_etalons,
)


# Соответствие нуклид → представитель цепочки для дочерних линий.
# МЭС решает одну A_k на эталон (цепочка в равновесии). Для отчёта
# дочерние нуклиды наследуют A родителя цепочки.
CHAIN_DAUGHTERS: Dict[str, List[str]] = {
    "Th-232": ["Ac-228", "Pb-212", "Tl-208", "Bi-212"],
    "Ra-226": ["Pb-214", "Bi-214"],
}


def _spec_energies(spec) -> np.ndarray:
    n = int(spec.n_channels)
    return np.array([spec.channel_to_energy(ch) for ch in range(n)], dtype=float)


def build_matrix_guard(sample, infos) -> dict:
    """F-448 matrix-compatibility guard.

    МЭС корректен ТОЛЬКО когда матрица/плотность эталона совпадает с
    пробой (LSRM §12.1: счёт ∝ объёмной активности Бк/л = Бк/кг·ρ с
    самопоглощением, зависящим от ρ/состава). Если плотности различаются,
    нормировка R на удельную активность (Бк/кг) даёт систематическую
    ошибку ≈ (ρ_проба/ρ_эталон), скорректированную самопоглощением.

    Источник (anti-hallucination): matrix_descriptor() из
    spec.extras['lsrm_sample_density_g_cm3'] / ['lsrm_material'].
    """
    s_name, s_rho = matrix_descriptor(sample)
    mismatches = []
    for i in infos:
        e_rho = i.density_g_cm3
        ratio = None
        flagged = False
        if s_rho and e_rho and e_rho > 0:
            ratio = round(s_rho / e_rho, 4)
            # >10% расхождение по плотности → нарушение требования МЭС
            flagged = abs(ratio - 1.0) > 0.10
        mismatches.append({
            "nuclide": i.nuclide,
            "etalon_matrix": i.matrix_name,
            "etalon_density_g_cm3": e_rho,
            "sample_matrix": s_name,
            "sample_density_g_cm3": s_rho,
            "density_ratio_sample_over_etalon": ratio,
            "mismatch": flagged,
        })
    any_mismatch = any(m["mismatch"] for m in mismatches)
    return {
        "sample_matrix": s_name,
        "sample_density_g_cm3": s_rho,
        "any_density_mismatch": any_mismatch,
        "per_etalon": mismatches,
        "rule": (
            "LSRM §12.1: МЭС требует одинаковую матрицу/плотность эталона "
            "и пробы. Расхождение ρ → систематическая ошибка A_МЭС "
            "≈ ρ_проба/ρ_эталон (скорр. самопоглощением). НЕ исправляется "
            "в Бк/кг-нормировке — нужны эталоны в той же матрице."
        ),
    }


def run_template(
    *,
    sample_path: Path,
    bg_path: Path,
) -> dict:
    """Построить R из эталонов и решить assay по пробе. Возвращает dict."""
    # ── Циркулярность ────────────────────────────────────────────────
    th_etalon = (REPO / ETALON_FILES["Th-232"]).resolve()
    if sample_path.resolve() == th_etalon:
        raise SystemExit(
            "ЦИРКУЛЯРНОСТЬ: путь пробы совпадает с Th-эталоном "
            f"420-17031 ({th_etalon}). Проба должна быть независимым "
            "файлом (для поверки — 420-7-17, A=1940 Бк/кг)."
        )

    R, infos, bg_meta = build_sensitivity_matrix_from_etalons(
        repo_root=REPO, bg_file=str(bg_path.relative_to(REPO))
        if bg_path.is_relative_to(REPO) else str(bg_path),
    )

    sample = read_lsrm_spe(str(sample_path))
    bg = read_lsrm_spe(str(bg_path))

    sample_counts = tuple(float(c) for c in np.asarray(sample.counts, dtype=float))
    sample_en = _spec_energies(sample)
    bg_counts = tuple(float(c) for c in np.asarray(bg.counts, dtype=float))
    bg_en = _spec_energies(bg)

    # template_assay вычитает окна фона внутри (bg-first, §12). Передаём
    # оси энергий пробы; фон — через BackgroundSpec + bg-энергии.
    bg_spec = BackgroundSpec(counts=bg_counts, t_live_s=float(bg.live_time))

    # template_assay использует единую energies_per_ch для пробы И фона.
    # Калибровки пробы/фона различаются (a0 sample=-8.29, bg=-5.04), но
    # в keV-domain окна инвариантны. Чтобы корректно посчитать окна фона
    # в ЕГО оси, пересоберём bg-rates через sum_in_windows на bg_en и
    # подставим как «эффективный» фон в сетке пробы. Проще: вызвать
    # template_assay с осью пробы, но фон уже учтён в R. Здесь bg
    # вычитается ещё раз для пробы — это корректно (R = (S-B)/A,
    # assay: (n-B) = R·A).
    result = template_assay(
        sample_counts=sample_counts,
        sample_t_live_s=float(sample.live_time),
        sensitivity=R,
        bg_spec=bg_spec,
        energies_per_ch=sample_en,
    )

    # F-448 matrix-compatibility guard (LSRM §12.1)
    matrix_guard = build_matrix_guard(sample, infos)

    rec = result.by_nuclide()
    nuclide_rows = []
    for nuc, (A, sA) in rec.items():
        nuclide_rows.append({
            "nuclide": nuc,
            "A_Bq_per_kg": round(float(A), 3),
            "sigma_Bq_per_kg": round(float(sA), 3),
            "is_chain_parent": True,
        })
        # дочерние нуклиды наследуют A родителя цепочки
        for d in CHAIN_DAUGHTERS.get(nuc, []):
            nuclide_rows.append({
                "nuclide": d,
                "A_Bq_per_kg": round(float(A), 3),
                "sigma_Bq_per_kg": round(float(sA), 3),
                "is_chain_parent": False,
                "inherits_from": nuc,
            })

    return {
        "method": "template_MES_LSRM_§12",
        "sample_path": str(sample_path),
        "bg": bg_meta,
        "windows_keV": [
            {"lo": lo, "hi": hi, "provenance": prov}
            for lo, hi, prov in MES_WINDOWS
        ],
        "R_matrix": {
            "nuclides": list(R.nuclides),
            "shape": list(R.R.shape),
            "values": [[round(float(x), 6) for x in row] for row in R.R],
        },
        "etalons": [
            {
                "nuclide": i.nuclide, "path": i.path,
                "A_certified_Bq_per_kg": i.A_certified_Bq_per_kg,
                "sigma_A_rel": i.sigma_A_rel, "t_live_s": i.t_live_s,
                "cert_raw": i.cert_raw,
            }
            for i in infos
        ],
        "chi2_per_dof": round(float(result.chi2_per_dof), 4),
        "n_windows": result.n_windows,
        "converged": result.converged,
        "notes": result.notes,
        "matrix_guard": matrix_guard,
        "nuclides": nuclide_rows,
    }


def load_efficiency_results(eff_json: Path) -> Dict[str, dict]:
    """Извлечь A_efficiency по нуклидам из cert_zcheck.json ИЛИ report.json.

    cert_zcheck.json: {"nuclides":[{"nuclide","A_Bq_per_kg","sigma_Bq_per_kg"}]}
    report.json:      {"identified_nuclides":[{"nuclide",
                       "specific_activity_Bq_per_kg",
                       "specific_activity_sigma_Bq_per_kg"}]}
    """
    d = json.loads(eff_json.read_text(encoding="utf-8"))
    out: Dict[str, dict] = {}
    if "nuclides" in d and d.get("nuclides"):
        for r in d["nuclides"]:
            nuc = r.get("nuclide")
            if nuc is None:
                continue
            out[nuc] = {
                "A_Bq_per_kg": r.get("A_Bq_per_kg"),
                "sigma_Bq_per_kg": r.get("sigma_Bq_per_kg"),
            }
    elif "identified_nuclides" in d:
        for r in d["identified_nuclides"]:
            nuc = r.get("nuclide")
            if nuc is None:
                continue
            out[nuc] = {
                "A_Bq_per_kg": r.get("specific_activity_Bq_per_kg"),
                "sigma_Bq_per_kg": r.get("specific_activity_sigma_Bq_per_kg"),
            }
    return out


def _delta_z(A: Optional[float], sA: Optional[float],
             A_ref: float, sA_ref: float):
    """Δ% и |z| против reference (паспорт ИЛИ другой метод)."""
    if A is None or A_ref is None or A_ref == 0:
        return None, None
    delta_pct = (A - A_ref) / A_ref * 100.0
    s = math.sqrt((sA or 0.0) ** 2 + (sA_ref or 0.0) ** 2)
    z = (A - A_ref) / s if s > 0 else None
    return round(delta_pct, 3), (round(abs(z), 3) if z is not None else None)


def build_compare(
    template_result: dict,
    eff_results: Dict[str, dict],
    *,
    cert_A: float,
    cert_rel: float,
    cert_name: str,
) -> dict:
    """Собрать таблицу A_eff vs A_МЭС vs паспорт с Δ%/|z|."""
    sig_cert = cert_A * cert_rel
    tmpl = {r["nuclide"]: r for r in template_result["nuclides"]}
    all_nuc = sorted(set(tmpl) | set(eff_results))

    rows = []
    for nuc in all_nuc:
        t = tmpl.get(nuc)
        e = eff_results.get(nuc)
        A_t = t["A_Bq_per_kg"] if t else None
        sA_t = t["sigma_Bq_per_kg"] if t else None
        A_e = e["A_Bq_per_kg"] if e else None
        sA_e = e["sigma_Bq_per_kg"] if e else None

        # Δ%/|z| против паспорта (применимо к нуклиду с тем же сертификатом —
        # для Th-232 chain паспорт = 1940; для прочих нет паспорта пробы).
        dt_cert_pct, dt_cert_z = _delta_z(A_t, sA_t, cert_A, sig_cert)
        de_cert_pct, de_cert_z = _delta_z(A_e, sA_e, cert_A, sig_cert)
        # МЭС vs efficiency (метод-к-методу)
        m2m_pct, m2m_z = _delta_z(A_t, sA_t, A_e, sA_e) if (A_e) else (None, None)

        rows.append({
            "nuclide": nuc,
            "A_efficiency_Bq_per_kg": A_e,
            "sigma_efficiency": sA_e,
            "A_template_MES_Bq_per_kg": A_t,
            "sigma_template": sA_t,
            "is_chain_parent": t.get("is_chain_parent") if t else None,
            "vs_cert": {
                "efficiency": {"delta_pct": de_cert_pct, "abs_z": de_cert_z},
                "template": {"delta_pct": dt_cert_pct, "abs_z": dt_cert_z},
            },
            "template_vs_efficiency": {"delta_pct": m2m_pct, "abs_z": m2m_z},
        })

    return {
        "certificate": {
            "name": cert_name, "A_Bq_per_kg": cert_A,
            "sigma_relative": cert_rel, "sigma_Bq_per_kg": round(sig_cert, 3),
        },
        "note": (
            "Δ%/|z| vs cert применимы только к нуклидам с известным "
            "паспортом пробы (Th-232 chain для Th-232 ОИСН-16). Для прочих "
            "нуклидов паспорт пробы отсутствует — сравнение метод-к-методу."
        ),
        "rows": rows,
    }


def _matrix_guard_md(guard: Optional[dict]) -> list:
    """F-448: markdown-блок matrix-guard (предупреждение о несовпадении ρ)."""
    if not guard:
        return []
    out = ["## Matrix-guard (совместимость матрицы эталон↔проба, LSRM §12.1)", ""]
    if guard.get("any_density_mismatch"):
        out.append("> **ВНИМАНИЕ: матрица/плотность пробы НЕ совпадает с эталонами.**")
        out.append("> Результат МЭС систематически смещён — см. правило ниже.")
    else:
        out.append("> Матрица/плотность эталонов и пробы совместимы (≤10%).")
    out.append("")
    out.append("| Эталон | Матрица эт. | ρ эт., г/см³ | Матрица пробы | ρ пробы | ρ_пр/ρ_эт | Несовпадение |")
    out.append("|---|---|---:|---|---:|---:|:---:|")
    for m in guard.get("per_etalon", []):
        flag = "ДА" if m.get("mismatch") else "—"
        out.append(
            f"| {m['nuclide']} | {m['etalon_matrix'] or '—'} "
            f"| {m['etalon_density_g_cm3'] if m['etalon_density_g_cm3'] is not None else '—'} "
            f"| {m['sample_matrix'] or '—'} "
            f"| {m['sample_density_g_cm3'] if m['sample_density_g_cm3'] is not None else '—'} "
            f"| {m['density_ratio_sample_over_etalon'] if m['density_ratio_sample_over_etalon'] is not None else '—'} "
            f"| {flag} |"
        )
    out.append("")
    out.append(f"> {guard.get('rule','')}")
    out.append("")
    return out


def compare_to_markdown(compare: dict, template_result: dict) -> str:
    cert = compare["certificate"]
    lines = []
    lines.append("# F-448 — Сравнение efficiency vs МЭС (template) vs паспорт")
    lines.append("")
    lines.append(f"**Проба:** `{template_result['sample_path']}`")
    lines.append(f"**Паспорт:** {cert['name']} = {cert['A_Bq_per_kg']:.0f} "
                 f"± {cert['sigma_relative']*100:.0f}% Бк/кг")
    lines.append(f"**χ²/dof МЭС:** {template_result['chi2_per_dof']}  "
                 f"(окон={template_result['n_windows']}, "
                 f"notes: {template_result['notes'] or '—'})")
    lines.append("")
    lines.extend(_matrix_guard_md(template_result.get("matrix_guard")))
    lines.append("## Эталоны (источник матрицы R)")
    lines.append("")
    lines.append("| Нуклид | Файл | A_серт, Бк/кг | σ% | t_live, с |")
    lines.append("|---|---|---:|---:|---:|")
    for et in template_result["etalons"]:
        fn = Path(et["path"]).name
        lines.append(
            f"| {et['nuclide']} | {fn} | {et['A_certified_Bq_per_kg']:.0f} | "
            f"{et['sigma_A_rel']*100:.0f} | {et['t_live_s']:.0f} |"
        )
    lines.append("")
    lines.append("## Окна МЭС (провенанс по линиям, data/nuclides.json)")
    lines.append("")
    lines.append("| Окно, кэВ | Линия / нуклид |")
    lines.append("|---|---|")
    for w in template_result["windows_keV"]:
        lines.append(f"| {w['lo']:.0f}–{w['hi']:.0f} | {w['provenance']} |")
    lines.append("")
    lines.append("## Сравнение активностей")
    lines.append("")
    lines.append(
        "| Нуклид | A_eff, Бк/кг | A_МЭС, Бк/кг | Δ%(МЭС vs eff) | "
        "|z|(МЭС vs eff) | Δ%(eff vs паспорт) | Δ%(МЭС vs паспорт) |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")

    def _fmt(x, suf=""):
        return f"{x:.1f}{suf}" if isinstance(x, (int, float)) else "—"

    for r in compare["rows"]:
        m2m = r["template_vs_efficiency"]
        lines.append(
            f"| {r['nuclide']} "
            f"| {_fmt(r['A_efficiency_Bq_per_kg'])} "
            f"| {_fmt(r['A_template_MES_Bq_per_kg'])} "
            f"| {_fmt(m2m['delta_pct'], '%') if m2m['delta_pct'] is not None else '—'} "
            f"| {_fmt(m2m['abs_z']) if m2m['abs_z'] is not None else '—'} "
            f"| {_fmt(r['vs_cert']['efficiency']['delta_pct'], '%') if r['vs_cert']['efficiency']['delta_pct'] is not None else '—'} "
            f"| {_fmt(r['vs_cert']['template']['delta_pct'], '%') if r['vs_cert']['template']['delta_pct'] is not None else '—'} |"
        )
    lines.append("")
    lines.append(f"> {compare['note']}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="F-448 МЭС/template assay + compare")
    ap.add_argument("--sample", type=str, default=str(REPO / SAMPLE_TH232_FILE))
    ap.add_argument("--bg", type=str, default=str(REPO / DEFAULT_BG_FILE))
    ap.add_argument("--mass", type=float, default=1.6)
    ap.add_argument("--efficiency-json", type=str, default=None)
    ap.add_argument("--cert-A", type=float, default=1940.0)
    ap.add_argument("--cert-rel", type=float, default=0.06)
    ap.add_argument("--cert-name", type=str,
                    default="Th-232 ОИСН-16 17-09-2007, m=1.6 кг")
    ap.add_argument("--out-dir", type=str, default=None)
    args = ap.parse_args()

    sample_path = Path(args.sample).resolve()
    bg_path = Path(args.bg).resolve()
    if not sample_path.exists():
        raise SystemExit(f"FAIL: проба не найдена: {sample_path}")
    if not bg_path.exists():
        raise SystemExit(f"FAIL: фон не найден: {bg_path}")

    print(f"[mes] sample: {sample_path}")
    print(f"[mes] bg:     {bg_path}")
    template_result = run_template(sample_path=sample_path, bg_path=bg_path)

    print(f"[mes] R.shape={template_result['R_matrix']['shape']} "
          f"χ²/dof={template_result['chi2_per_dof']} "
          f"converged={template_result['converged']}")
    print("[mes] A_k(МЭС):")
    for r in template_result["nuclides"]:
        if r.get("is_chain_parent"):
            print(f"   {r['nuclide']:8s} = {r['A_Bq_per_kg']:.1f} "
                  f"± {r['sigma_Bq_per_kg']:.1f} Бк/кг")

    mg = template_result.get("matrix_guard") or {}
    if mg.get("any_density_mismatch"):
        print("[mes] *** MATRIX-GUARD: матрица/плотность пробы НЕ совпадает "
              "с эталонами — A_МЭС систематически смещён! ***")
        for m in mg.get("per_etalon", []):
            if m.get("mismatch"):
                print(f"   {m['nuclide']:8s}: эталон {m['etalon_matrix']} "
                      f"ρ={m['etalon_density_g_cm3']} vs проба "
                      f"{m['sample_matrix']} ρ={m['sample_density_g_cm3']} "
                      f"(ρ_пр/ρ_эт={m['density_ratio_sample_over_etalon']})")

    # ── out-dir ──────────────────────────────────────────────────────
    if args.out_dir:
        out_dir = Path(args.out_dir).resolve()
    elif args.efficiency_json:
        out_dir = Path(args.efficiency_json).resolve().parent
    else:
        ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        out_dir = (REPO / "demo_reports" / f"{ts}_template_assay").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    eff_results: Dict[str, dict] = {}
    if args.efficiency_json:
        eff_path = Path(args.efficiency_json).resolve()
        if eff_path.exists():
            eff_results = load_efficiency_results(eff_path)
            print(f"[mes] efficiency results loaded: {len(eff_results)} nuclides "
                  f"from {eff_path.name}")
        else:
            print(f"[mes] WARN: efficiency-json не найден: {eff_path}")

    compare = build_compare(
        template_result, eff_results,
        cert_A=args.cert_A, cert_rel=args.cert_rel, cert_name=args.cert_name,
    )

    # ── write artifacts ──────────────────────────────────────────────
    out_json = out_dir / "template_compare.json"
    payload = {
        "template_assay": template_result,
        "compare": compare,
        "mass_kg": args.mass,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    out_json.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    out_md = out_dir / "template_compare.md"
    out_md.write_text(
        compare_to_markdown(compare, template_result), encoding="utf-8",
    )

    print()
    print("=" * 72)
    print("F-448 АРТЕФАКТЫ (МЭС + compare):")
    print(f"  • JSON: {out_json}")
    print(f"  • MD:   {out_md}")
    print("=" * 72)
    print(compare_to_markdown(compare, template_result))
    return 0


if __name__ == "__main__":
    sys.exit(main())