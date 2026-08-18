"""cert_zcheck — z-test + residual % confirmed nuclides vs certificate.

Постпроцесс после run_skill.py: читает `<bundle>/sample_v2/<stem>_report.json`,
для каждого confirmed nuclide считает:

    Δ_% = (A - A_cert) / A_cert * 100
    z   = (A - A_cert) / sqrt(σ_A² + σ_cert²)

Phase 1 exit-criteria (CLAUDE.md → Phase 1 exit-criteria):
    PASS = |Δ_%| ≤ 15  AND  |z| ≤ 2

Пишет таблицу в stdout + JSON `<bundle>/cert_zcheck.json`.

Usage:
    PYTHONIOENCODING=utf-8 python scripts/cert_zcheck.py <bundle-dir> \
        --cert-A 1940 --cert-rel 0.06 [--cert-name "Th-232 ОИСН-16"]

Используется автоматически из `scripts/run_plan_a.py`, если переданы
defaults / env vars `GAMMA_CERT_A`, `GAMMA_CERT_REL`.

Зафиксировано оператором 2026-06-10: «почему на фоне Z не вычислено? Переделай.
Зафиксируй». Process-bug: Plan A summary до v1.2.7 не вычислял |z|, только Δ%.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def find_report_json(bundle: Path) -> Path:
    v2 = sorted((bundle / "sample_v2").glob("*_report.json"))
    if v2:
        return v2[0]
    prod = sorted((bundle / "sample").glob("*_report.json"))
    if prod:
        return prod[0]
    raise SystemExit(f"FAIL: no *_report.json found in {bundle}/sample_v2/ or {bundle}/sample/")


def compute(report_json: Path, A_cert: float, rel_cert: float, cert_name: str):
    d = json.loads(report_json.read_text(encoding="utf-8"))
    sig_cert = A_cert * rel_cert
    out = {
        "certificate": {
            "name": cert_name,
            "A_Bq_per_kg": A_cert,
            "sigma_relative": rel_cert,
            "sigma_Bq_per_kg": round(sig_cert, 3),
        },
        "phase1_criteria": {
            "residual_pct_max": 15.0,
            "abs_z_max": 2.0,
            "rule": "PASS if |Δ_%| ≤ 15 AND |z| ≤ 2",
        },
        "nuclides": [],
    }
    for n in d.get("identified_nuclides", []):
        A = n.get("specific_activity_Bq_per_kg")
        sA = n.get("specific_activity_sigma_Bq_per_kg")
        nu = n.get("nuclide")
        row = {"nuclide": nu, "tier": n.get("tier"), "characteristic_line_keV": n.get("characteristic_line_keV")}
        if A is None or sA is None:
            row.update({
                "A_Bq_per_kg": None, "sigma_Bq_per_kg": None,
                "delta_pct": None, "abs_z": None,
                "verdict": "n/a",
                "passes_residual": None, "passes_z": None,
            })
        else:
            delta_pct = (A - A_cert) / A_cert * 100.0
            sig_comb = math.sqrt(sA * sA + sig_cert * sig_cert)
            z = (A - A_cert) / sig_comb
            pass_r = abs(delta_pct) <= 15.0
            pass_z = abs(z) <= 2.0
            if pass_r and pass_z:
                verdict = "PASS"
            elif not pass_z:
                verdict = "FAIL |z|>2"
            else:
                verdict = "FAIL Δ>15%"
            row.update({
                "A_Bq_per_kg": round(A, 3),
                "sigma_Bq_per_kg": round(sA, 3),
                "delta_pct": round(delta_pct, 3),
                "sigma_combined_Bq_per_kg": round(sig_comb, 3),
                "z": round(z, 4),
                "abs_z": round(abs(z), 4),
                "passes_residual": pass_r,
                "passes_z": pass_z,
                "verdict": verdict,
            })
        out["nuclides"].append(row)
    confirmed_count = sum(1 for r in out["nuclides"] if r.get("verdict") not in ("n/a", "upper-limit"))
    pass_count = sum(1 for r in out["nuclides"] if r.get("verdict") == "PASS")
    out["summary"] = {
        "n_confirmed": confirmed_count,
        "n_passed": pass_count,
        "n_failed": confirmed_count - pass_count,
        "all_pass": confirmed_count > 0 and pass_count == confirmed_count,
    }
    return out


def print_table(out: dict):
    cert = out["certificate"]
    print(f"Cert: {cert['name']} = {cert['A_Bq_per_kg']:.1f} ± {cert['sigma_relative']*100:.1f}% Бк/кг "
          f"(σ_cert = {cert['sigma_Bq_per_kg']:.1f} Бк/кг)")
    print(f"Phase 1 criteria: |Δ_%| ≤ 15 AND |z| ≤ 2 на combined stat+syst")
    print()
    print(f"{'Nuclide':<10}{'A,Бк/кг':>10}{'σ_meas':>10}{'Δ_%':>10}{'σ_comb':>10}{'|z|':>8}  Verdict")
    print("-" * 72)
    for r in out["nuclides"]:
        if r["A_Bq_per_kg"] is None:
            print(f"{r['nuclide']:<10}{'<DL':>10}{'-':>10}{'-':>10}{'-':>10}{'-':>8}  {r['verdict']}")
        else:
            print(f"{r['nuclide']:<10}{r['A_Bq_per_kg']:>10.1f}{r['sigma_Bq_per_kg']:>10.1f}"
                  f"{r['delta_pct']:>+10.2f}{r['sigma_combined_Bq_per_kg']:>10.1f}{r['abs_z']:>8.2f}  {r['verdict']}")
    s = out["summary"]
    print()
    print(f"Summary: {s['n_passed']}/{s['n_confirmed']} PASS"
          + (" — ALL PASS Phase 1 exit-criteria" if s["all_pass"] else " — FAIL"))


def main():
    ap = argparse.ArgumentParser(description="z-test + Δ% vs certificate, Phase 1 exit-criteria")
    ap.add_argument("bundle", type=Path, help="Path to run_skill.py output bundle")
    ap.add_argument("--cert-A", type=float, required=True, help="Certified activity, Бк/кг")
    ap.add_argument("--cert-rel", type=float, required=True, help="Certified relative sigma (e.g. 0.06 for ±6%%)")
    ap.add_argument("--cert-name", type=str, default="certificate", help="Cert description")
    ap.add_argument("--out", type=Path, default=None, help="Output JSON path (default: <bundle>/cert_zcheck.json)")
    args = ap.parse_args()

    bundle = args.bundle.resolve()
    if not bundle.is_dir():
        raise SystemExit(f"FAIL: bundle not a directory: {bundle}")
    report = find_report_json(bundle)
    out = compute(report, args.cert_A, args.cert_rel, args.cert_name)
    out["source_report"] = str(report)
    out_path = args.out or (bundle / "cert_zcheck.json")
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print_table(out)
    print()
    print(f"JSON: {out_path}")
    return 0 if out["summary"]["all_pass"] else 3


if __name__ == "__main__":
    sys.exit(main())