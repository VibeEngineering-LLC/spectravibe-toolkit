"""diag_vs_lsrm_th232_marinelli.py

Комплексная сверка результата нашего pipeline (V2 JSON-отчёт run_plan_a.py)
с эталоном LSRM Гамма-1С для Th232_420-7-17_Маринелли_0cm.

Эталон: references/lsrm_ground_truth/Th232_420-7-17_Marinelli_0cm/
  peaks_table_lsrm.json
  fwhm_calibration_lsrm.json
  peak_zones_lsrm.json
  activity_lsrm.json

Сверка по 5 осям (README протокол §43-66 эталона):
  1. Энергии пиков (Δ E)
  2. ПШПВ (FWHM точечно + кривая)
  3. Площади (counts)
  4. Активности (Бк/кг)
  5. Параметры формы (bg_polynomial_degree, step_value sign)

Usage:
  PYTHONIOENCODING=utf-8 python scripts/diag/diag_vs_lsrm_th232_marinelli.py \
      --report demo_reports/<ts>_<label>/sample_v2/Th232_*_report.json

При наличии --report по умолчанию сохраняет страницу сравнения рядом с отчётом:
  <report_dir>/comparison_vs_lsrm.html  (для оператора, открыть в браузере)
  <report_dir>/comparison_vs_lsrm.md    (raw markdown, для git/diff)
Подавить файловый вывод можно флагом --stdout-only.

Если --report не задан — берёт самый свежий из demo_reports/*th232*marinelli*/sample_v2/*_report.json.
"""

import argparse
import glob
import html as html_lib
import os
import re
import sys
import json
from pathlib import Path


REPO = Path(r"<WORKDIR>\gamma-spectrum-analysis")
LSRM_DIR = REPO / "references" / "lsrm_ground_truth" / "Th232_420-7-17_Marinelli_0cm"
SAMPLE_MASS_KG = 1.6
CERT_BQ_PER_KG = 1940.0
CERT_SIGMA_REL = 0.06

ENERGY_MATCH_TOL_KEV = 5.0
ENERGY_PASS_TOL_KEV = 0.5
FWHM_PASS_REL = 0.10
AREA_PASS_SIGMA = 2.0


def find_latest_report():
    pattern = str(REPO / "demo_reports" / "*th232*marinelli*" / "sample_v2" / "*report.json")
    cands = sorted(glob.glob(pattern), key=lambda p: os.path.getmtime(p))
    if not cands:
        raise SystemExit(f"[FATAL] нет отчётов по шаблону {pattern}")
    return cands[-1]


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def fwhm_model_value(coeffs, E):
    a, b, c = (coeffs + [0.0, 0.0, 0.0])[:3]
    sq = a + b * E + c * E * E
    return (sq ** 0.5) if sq > 0 else float("nan")


def classify_delta(delta_abs, tol_pass, tol_marginal_mult=2.0):
    if delta_abs <= tol_pass:
        return "PASS"
    if delta_abs <= tol_pass * tol_marginal_mult:
        return "MARGINAL"
    return "FAIL"


def match_by_lib_energy(our_feps, lsrm_peaks):
    matches = []
    our_by_lib = {}
    for fep in our_feps:
        libE = fep.get("library_E_keV")
        if libE is None:
            continue
        our_by_lib.setdefault(round(libE, 2), []).append(fep)
    for lp in lsrm_peaks:
        libE = lp["library_energy_keV"]
        best = None
        best_dE = 1e9
        for k, feps in our_by_lib.items():
            d = abs(k - libE)
            if d < best_dE and d <= ENERGY_MATCH_TOL_KEV:
                best_dE = d
                best = feps[0]
        matches.append((lp, best))
    return matches


# ------------------ markdown → html минималистичный конвертер ------------------
_VERDICT_CLASS = {
    "PASS": "v-pass", "MARGINAL": "v-marg", "FAIL": "v-fail", "MISS": "v-miss",
}


def _md_inline(s: str) -> str:
    out = html_lib.escape(s, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", out)
    # Подсветка verdict-токенов в строках
    for v, cls in _VERDICT_CLASS.items():
        out = re.sub(rf"<strong>{v}</strong>", f'<strong class="{cls}">{v}</strong>', out)
        out = re.sub(rf"(?<![A-Za-z]){v}(?![A-Za-z])", f'<span class="{cls}">{v}</span>', out)
    return out


def md_to_html(md: str, *, title: str) -> str:
    body_parts = []
    lines = md.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line:
            i += 1
            continue
        if line.startswith("---"):
            body_parts.append("<hr/>")
            i += 1
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            level = len(m.group(1))
            body_parts.append(f"<h{level}>{_md_inline(m.group(2))}</h{level}>")
            i += 1
            continue
        # table: header + separator + rows
        if line.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s:\-|]+\|$", lines[i + 1].rstrip()):
            header = [c.strip() for c in line.strip("|").split("|")]
            i += 2
            rows = []
            while i < len(lines) and lines[i].rstrip().startswith("|"):
                rows.append([c.strip() for c in lines[i].rstrip().strip("|").split("|")])
                i += 1
            thead = "".join(f"<th>{_md_inline(h)}</th>" for h in header)
            tbody = ""
            for r in rows:
                # row-level verdict подсветка (по последнему столбцу)
                last_token = r[-1] if r else ""
                row_class = ""
                for v, cls in _VERDICT_CLASS.items():
                    if f"**{v}**" in last_token or f">{v}<" in last_token or f" {v} " in (" " + last_token + " "):
                        row_class = f' class="row-{cls}"'
                        break
                cells = "".join(f"<td>{_md_inline(c)}</td>" for c in r)
                tbody += f"<tr{row_class}>{cells}</tr>"
            body_parts.append(
                f'<table class="diag"><thead><tr>{thead}</tr></thead><tbody>{tbody}</tbody></table>'
            )
            continue
        # list bullet
        if line.startswith("- "):
            items = []
            while i < len(lines) and lines[i].rstrip().startswith("- "):
                items.append(f"<li>{_md_inline(lines[i].rstrip()[2:])}</li>")
                i += 1
            body_parts.append("<ul>" + "".join(items) + "</ul>")
            continue
        # paragraph
        body_parts.append(f"<p>{_md_inline(line)}</p>")
        i += 1
    body = "\n".join(body_parts)
    return f"""<!DOCTYPE html>
<html lang="ru"><head><meta charset="utf-8"><title>{html_lib.escape(title)}</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ font-family: -apple-system, Segoe UI, Roboto, sans-serif; max-width: 1200px;
       margin: 1.5rem auto; padding: 0 1rem; line-height: 1.45; }}
h1 {{ border-bottom: 2px solid #888; padding-bottom: .3em; }}
h2 {{ border-bottom: 1px solid #888; padding-bottom: .2em; margin-top: 2em; }}
h3 {{ margin-top: 1.5em; }}
table.diag {{ border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.92em; }}
table.diag th, table.diag td {{ border: 1px solid #888; padding: 4px 8px; text-align: right; }}
table.diag th {{ background: rgba(127,127,127,.15); }}
table.diag td:first-child, table.diag th:first-child {{ text-align: left; }}
table.diag td:last-child, table.diag th:last-child {{ text-align: left; white-space: nowrap; }}
code {{ background: rgba(127,127,127,.18); padding: 1px 4px; border-radius: 3px;
       font-size: 0.9em; font-family: Consolas, monospace; }}
.v-pass {{ color: #1b8a3a; font-weight: bold; }}
.v-marg {{ color: #b07b00; font-weight: bold; }}
.v-fail {{ color: #c1342f; font-weight: bold; }}
.v-miss {{ color: #6a6a6a; font-style: italic; font-weight: bold; }}
tr.row-v-pass td {{ background: rgba(27,138,58,.07); }}
tr.row-v-marg td {{ background: rgba(176,123,0,.10); }}
tr.row-v-fail td {{ background: rgba(193,52,47,.10); }}
tr.row-v-miss td {{ background: rgba(127,127,127,.06); color: #555; }}
hr {{ border: none; border-top: 1px dashed #888; margin: 2em 0; }}
ul {{ margin: 0.4em 0 0.8em 1.2em; }}
</style></head><body>
<h1>{html_lib.escape(title)}</h1>
{body}
</body></html>"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--report", default=None, help="Путь к V2 JSON-отчёту run_plan_a.py")
    p.add_argument("--out-dir", default=None,
                   help="Куда сохранить comparison_vs_lsrm.{html,md} (по умолчанию рядом с отчётом)")
    p.add_argument("--stdout-only", action="store_true",
                   help="Не писать файлы, только stdout markdown")
    args = p.parse_args()

    report_path = args.report or find_latest_report()
    out_lines = []

    def w(s=""):
        out_lines.append(s)

    w(f"## Источник нашего отчёта\n`{report_path}`\n")
    rep = load(report_path)

    lsrm_peaks = load(LSRM_DIR / "peaks_table_lsrm.json")["peaks"]
    lsrm_fwhm = load(LSRM_DIR / "fwhm_calibration_lsrm.json")
    lsrm_zones = load(LSRM_DIR / "peak_zones_lsrm.json")["zones"]
    lsrm_act = load(LSRM_DIR / "activity_lsrm.json")

    feps = rep.get("primary_feps") or []
    ids = rep.get("identified_nuclides") or []
    cal = rep.get("calibration") or {}
    fwhm_cal = cal.get("fwhm_cal") or {}
    fwhm_src = fwhm_cal.get("source", "?")
    fwhm_coef = fwhm_cal.get("coefficients") or [0.0, 0.0, 0.0]

    w(f"## Сводка нашего отчёта")
    w(f"- primary_feps: {len(feps)} линий")
    w(f"- identified_nuclides: {len(ids)} нуклидов ({', '.join(n['nuclide'] for n in ids)})")
    w(f"- fwhm_cal.source: `{fwhm_src}` (LSRM = polynomial deg=4, χ²=2.0332)")
    w(f"- fwhm_cal.coefficients: {fwhm_coef}")
    w()

    w("---")
    w("## 1. Энергии и ПШПВ пиков (matched по library_E)")
    w()
    w("| LSRM lib_E | LSRM E_max | наш E_max | Δ E | LSRM FWHM | наш FWHM | Δ FWHM | LSRM area | наш area | Δ area | verdict |")
    w("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|")

    matches = match_by_lib_energy(feps, lsrm_peaks)
    counts = {"PASS": 0, "MARGINAL": 0, "FAIL": 0, "MISS": 0}
    fail_lines = []

    for lp, our in matches:
        libE = lp["library_energy_keV"]
        lsrm_E = lp["energy_keV"]
        lsrm_F = lp["fwhm_keV"]
        lsrm_A = lp["area_counts"]
        lsrm_sA = lp["area_uncertainty_counts"]

        if our is None:
            w(f"| {libE:.2f} | {lsrm_E:.2f} | — | — | {lsrm_F:.2f} | — | — | {lsrm_A:.0f} | — | — | **MISS** |")
            counts["MISS"] += 1
            fail_lines.append(f"MISS: библиотечная линия {libE} кэВ — нет matched FEP в нашем primary_feps")
            continue

        our_E = our.get("peak_E_keV")
        our_F = our.get("fwhm_keV")
        our_A = our.get("peak_area_counts")
        our_sA = our.get("peak_area_sigma")

        dE = our_E - lsrm_E if our_E is not None else None
        dF = our_F - lsrm_F if our_F is not None else None
        dA = our_A - lsrm_A if our_A is not None else None
        tol_E = max(lp.get("energy_uncertainty_keV", 0.5), ENERGY_PASS_TOL_KEV)
        verdict_E = classify_delta(abs(dE), tol_E) if dE is not None else "FAIL"
        tol_F = max(lp.get("fwhm_uncertainty_meas_eV", 0) / 1000.0, 2.0)
        verdict_F = classify_delta(abs(dF), tol_F) if dF is not None else "FAIL"
        if dA is not None and (our_sA or lsrm_sA):
            sig_combo = max(our_sA or 0.0, lsrm_sA or 0.0, 1.0)
            n_sig = abs(dA) / sig_combo
            if n_sig <= AREA_PASS_SIGMA:
                verdict_A = "PASS"
            elif n_sig <= 3.0:
                verdict_A = "MARGINAL"
            else:
                verdict_A = "FAIL"
        else:
            verdict_A = "FAIL"

        joint = "FAIL" if "FAIL" in (verdict_E, verdict_F, verdict_A) else (
            "MARGINAL" if "MARGINAL" in (verdict_E, verdict_F, verdict_A) else "PASS"
        )
        counts[joint] += 1
        if joint == "FAIL":
            fail_lines.append(
                f"FAIL E={libE}: E={verdict_E}(Δ{dE:+.2f}) F={verdict_F}(Δ{dF:+.2f}) A={verdict_A}(Δ{dA:+.0f})"
            )

        w(f"| {libE:.2f} | {lsrm_E:.2f} | {our_E:.2f} | {dE:+.2f} | "
          f"{lsrm_F:.2f} | {our_F:.2f} | {dF:+.2f} | "
          f"{lsrm_A:.0f} | {our_A:.0f} | {dA:+.0f} | "
          f"E={verdict_E} F={verdict_F} A={verdict_A} → **{joint}** |")

    w()
    w(f"**Итог matched:** PASS={counts['PASS']}, MARGINAL={counts['MARGINAL']}, FAIL={counts['FAIL']}, MISS={counts['MISS']} (всего {len(lsrm_peaks)} LSRM-линий)")
    w()

    w("---")
    w("## 2. Кривая FWHM(E) на 14 LSRM-якорях")
    w()
    w(f"Наша модель: source=`{fwhm_src}`, coef={fwhm_coef}")
    w()
    w("| E, keV | LSRM FWHM (изм) | LSRM FWHM (калибр) | наш FWHM(E) | Δ нашего | Δ % | tol(σ_meas) | verdict |")
    w("|---:|---:|---:|---:|---:|---:|---:|:---|")
    fwhm_pass = fwhm_marg = fwhm_fail = 0
    for a in lsrm_fwhm["anchors"]:
        E = a["energy_keV"]
        Fm = a["fwhm_keV_measured"]
        Fc = a["fwhm_keV_calibration"]
        our_F = fwhm_model_value(fwhm_coef, E)
        dF = our_F - Fm
        dF_pct = (dF / Fm) * 100 if Fm else 0
        tol = max(a.get("fwhm_uncertainty_meas_eV", 0) / 1000.0, 2.0)
        v = classify_delta(abs(dF), tol)
        if v == "PASS":
            fwhm_pass += 1
        elif v == "MARGINAL":
            fwhm_marg += 1
        else:
            fwhm_fail += 1
        w(f"| {E:.2f} | {Fm:.2f} | {Fc:.2f} | {our_F:.2f} | {dF:+.2f} | {dF_pct:+.1f}% | {tol:.2f} | **{v}** |")
    w()
    w(f"**FWHM-кривая итог:** PASS={fwhm_pass}, MARGINAL={fwhm_marg}, FAIL={fwhm_fail} (из 14 якорей)")
    w()
    if fwhm_pass < 7:
        w(f"  ⚠ Меньше половины якорей в PASS — модель FWHM не соответствует LSRM-эталону.")
        w(f"  Источник нашей модели: `{fwhm_src}` — если это `default_*`, bootstrap не активировался.")
        w()

    w("---")
    w("## 3. Параметры зон (bg_poly_deg, step_value, χ²)")
    w()
    w("| LSRM E | LSRM bg_deg | LSRM step | LSRM χ² | наш bg_deg | наш step (erfc) | замечания |")
    w("|---:|---:|---:|---:|:---|:---|:---|")
    for z in lsrm_zones:
        comments = []
        if z["bg_polynomial_degree"] != 1:
            comments.append("наш всегда deg=1")
        if z["step_value"] < 0:
            comments.append("LSRM step<0; наш erfc всегда нисходящий")
        w(f"| {z['library_energy_keV']:.2f} | {z['bg_polynomial_degree']} | "
          f"{z['step_value']:+.5f} | {z['chi2']:.3f} | 1 (фикс.) | erfc-step (нисходящий) | "
          f"{'; '.join(comments) if comments else 'ok-shape'} |")
    w()

    w("---")
    w("## 4. Активности по нуклидам (Бк/кг)")
    w()
    w(f"Сертификат Th-232 ОИСН-16 17-09-2007: A_specific = {CERT_BQ_PER_KG:.1f} ± {CERT_BQ_PER_KG*CERT_SIGMA_REL:.1f} Бк/кг (m={SAMPLE_MASS_KG} кг)")
    w()
    w("| нуклид | A_наш Бк/кг | σ_наш | LSRM 2614/911/238/583 → A/m=1.6 | дельта от cert (%) | verdict |")
    w("|:---|---:|---:|:---|---:|:---|")
    lsrm_bq_per_kg_avg = {}
    for line in lsrm_act["lines"]:
        if not isinstance(line.get("activity"), (int, float)):
            continue
        e = line["library_energy_keV"]
        lsrm_bq_per_kg_avg[e] = line["activity"] / SAMPLE_MASS_KG
    for nuc in ids:
        if nuc.get("specific_activity_Bq_per_kg") is None:
            continue
        A = nuc["specific_activity_Bq_per_kg"]
        sA = nuc.get("specific_activity_sigma_Bq_per_kg") or 0
        nuclide = nuc["nuclide"]
        rel_cert = (A - CERT_BQ_PER_KG) / CERT_BQ_PER_KG * 100
        v = "PASS" if abs(rel_cert) <= 15 else ("MARGINAL" if abs(rel_cert) <= 30 else "FAIL")
        lsrm_str = "; ".join(f"{e:.0f}→{bq:.0f}" for e, bq in sorted(lsrm_bq_per_kg_avg.items()) if e in (238.632, 583.187, 911.204, 2614.511))
        w(f"| {nuclide} | {A:.1f} | {sA:.1f} | {lsrm_str} | {rel_cert:+.1f}% | **{v}** |")
    w()

    w("---")
    w("## 5. Общий вердикт и приоритеты исправлений")
    w()
    total = sum(counts.values())
    pct_ok = (counts["PASS"] / total * 100) if total else 0
    w(f"- matched-сводка: **{counts['PASS']}/{total} PASS ({pct_ok:.0f}%)**")
    w(f"- FWHM-кривая: **{fwhm_pass}/14 PASS** (источник модели: `{fwhm_src}`)")
    if fail_lines:
        w()
        w("### FAIL-инциденты (для разбора)")
        for fl in fail_lines:
            w(f"- {fl}")
    w()
    w("### Известные систематики (НЕ требуют немедленного фикса по решению F-160)")
    w("- LSRM использует асимметричный Гаусс + адаптивный bg-deg (1/2/3) → даёт лучшую χ²/dof на низких энергиях.")
    w("- Знак step_value: LSRM ±, наш erfc всегда нисходящий (на 580 кэВ это ~3% эффект на counts).")
    w("- F-160 (2026-06-20): 3-параметрическая NNLS-аппроксимация LSRM-полинома 4-й степени имеет систематику ±5-7 keV на якорях; это структурное ограничение API `build_fwhm_model` ((a,b,c)-tuple).")

    md_text = "\n".join(out_lines) + "\n"
    print(md_text)

    if not args.stdout_only:
        out_dir = Path(args.out_dir) if args.out_dir else Path(report_path).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        md_path = out_dir / "comparison_vs_lsrm.md"
        html_path = out_dir / "comparison_vs_lsrm.html"
        md_path.write_text(md_text, encoding="utf-8")
        title = "Сравнение с эталоном LSRM Гамма-1С — Th-232 Marinelli 0 cm"
        html_path.write_text(md_to_html(md_text, title=title), encoding="utf-8")
        print(f"[diag] Сохранено: {md_path}", file=sys.stderr)
        print(f"[diag] Сохранено: {html_path}", file=sys.stderr)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
