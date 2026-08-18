"""Skill vs canon audit — point-by-point check that canonical pipeline
implements every requirement of the gamma-spectrum-analysis skill.

Usage:
    PYTHONIOENCODING=utf-8 python scripts/audit/skill_vs_canon_audit.py
"""
from __future__ import annotations
import sys, json, re, pathlib
from dataclasses import dataclass
from fnmatch import fnmatch

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = pathlib.Path(r"<WORKDIR>\gamma-spectrum-analysis")
CANON = ROOT / "scripts" / "gamma"

@dataclass
class Check:
    id: str
    step: str
    requirement: str
    patterns: list
    note: str = ""

CHECKS = [
    Check("WP1", "principles", "no smoothing on displayed data points",
          [(r"no[_-]smooth|raw.*counts|smoothing.*off|do not smooth", "reporting/*.py")]),
    Check("WP2", "principles", "polynomial degree <= 4 (energy/FWHM/efficiency)",
          [(r"deg(ree)?\s*<=?\s*4|max[_-]?deg(ree)?\s*=\s*4|polyfit.*,\s*4\b", "calibration/*.py")]),
    Check("WP3", "principles", "deconv only for already-identified nuclides",
          [(r"identified.*only|confirmed.*nuclide|expected[_-]?components", "peaks/deconvolve.py")]),
    Check("WP4", "principles", "library-directed peak search (candidate-list-only)",
          [(r"library[_-]directed|candidate[_-]?list|targeted.*search", "identification/*.py")]),

    Check("S1a", "1", "live & real time both required for cps",
          [(r"live_time|t_live", "io/*.py"),
           (r"real_time|t_real", "io/*.py")]),
    Check("S1b", "1", "multi-format parsers (.spe/.chn/.n42/.mca/.txt/.csv/.spm)",
          [(r"\.spe|\.chn|\.n42|\.mca|\.spm", "io/*.py")]),
    Check("S1c", "1", "filename token parsing",
          [(r"filename.*token|parse[_-]?filename|filename[_-]hint", "io/*.py")]),
    Check("S1d", "1", "BG-spectrum extraction from multi-block (.spm/.n42)",
          [(r"background.*section|measurementClass|spm.*section", "io/*.py")]),
    Check("S1e", "1", "BG calibration independent from sample",
          [(r"bg.*independent.*cal|background.*own.*cal|recal.*background", "calibration/*.py")]),
    Check("S1f", "1", "display: cps log no-smooth title-meta BG-overlay",
          [(r"log[_-]?scale|yaxis.*log|background.*trace|overlay.*background", "reporting/*.py")]),

    Check("S2a", "2", "Pb K-XRF 73-87 keV check",
          [(r"73\.0|74\.97|84\.4|84\.9|87\.3|Pb.*K[-_]?X", "identification/*.py")]),
    Check("S2b", "2", "natural BG nuclides K40/Tl208/Ac228/Bi214/Pb212",
          [(r"K[-_]?40.*1460|Tl[-_]?208.*2614|Ac[-_]?228|Bi[-_]?214|Pb[-_]?212", "**/*.py")]),
    Check("S2c", "2", "continuum-level diagnostic",
          [(r"continuum.*level|bg.*cps.*total|baseline.*cps", "calibration/*.py")]),
    Check("S2d", "2", "environment classification at top of report",
          [(r"low[_-]?background|shielded|natural[_-]?background|environment.*classif", "**/*.py")]),

    Check("S3a", "3", "Mariscotti 2nd-derivative search in CHANNEL space",
          [(r"mariscotti|second[_-]?derivative|2nd[_-]?deriv", "peaks/*.py")]),
    Check("S3b", "3", "Currie L_C threshold with k_alpha 1.645",
          [(r"1\.645|k[_-]?alpha|Currie|L[_-]?C\b", "peaks/*.py")]),
    Check("S3c", "3", "rough FWHM 5-20 ch HPGe / 30-100 ch scint",
          [(r"fwhm.*ch.*\d+|rough.*fwhm|expected.*fwhm.*ch", "peaks/*.py")]),

    Check("S4a", "4", "relative resolution check",
          [(r"relative.*resolution|R_pct|resolution.*662", "calibration/detector_type.py")]),
    Check("S4b", "4", "intrinsic-signature catalog",
          [(r"138La|La[-_]?138|227Ac|intrinsic", "calibration/detector_type.py")]),

    Check("S5a", "5", "skip stored E-cal if residuals < 0.3*FWHM",
          [(r"0\.3.*fwhm|stored.*residual|residual.*0\.3", "calibration/stored_check.py")]),
    Check("S5b", "5", "invariant anchor catalog (Co60 K40 Tl208 Cs137 ann)",
          [(r"1460\.82|2614\.51|661\.66|1173|1332|annihilation|511", "calibration/*.py")]),
    Check("S5c", "5", "filename-prior anchor priority",
          [(r"filename.*prior|filename.*hint.*anchor", "calibration/*.py")]),
    Check("S5d", "5", "E-cal poly deg<=4, piecewise above",
          [(r"piecewise|max[_-]?deg.*4|degree.*<=?.*4", "calibration/energy_fit.py")]),

    Check("S6a", "6", "skip stored FWHM if within 5%",
          [(r"5\s*%|0\.05.*fwhm|stored.*fwhm.*ok", "calibration/fwhm_provider.py")]),
    Check("S6b", "6", "HPGe FWHM^2=a+bE+cE^2",
          [(r"fwhm\*\*2|fwhm_squared|a.*\+.*b.*E.*\+.*c.*E", "calibration/fwhm_fit.py")]),
    Check("S6c", "6", "scint FWHM=k*sqrt(E+alpha*E^2) non-proportionality",
          [(r"non[_-]?proportional|sqrt.*alpha|k.*sqrt.*E", "calibration/fwhm_fit.py")]),
    Check("S6d", "6", "CdZnTe left-tail params",
          [(r"czt|CdZnTe|left[_-]?tail", "calibration/*.py")]),

    Check("S7a", "7", "candidate list build",
          [(r"candidate[_-]?list|build[_-]?candidates", "identification/*.py")]),
    Check("S7b", "7", "characteristic-line approach",
          [(r"characteristic[_-]?line|lowest[_-]?MDA", "identification/*.py")]),
    Check("S7c", "7", "+/- 2sigma_E energy match",
          [(r"2\s*\*?\s*sigma[_-]?E|2sigma|tolerance.*sigma_E", "identification/*.py")]),
    Check("S7d", "7", "FWHM match within 10%",
          [(r"10\s*%|0\.10.*fwhm|fwhm.*match.*0\.1", "identification/*.py")]),
    Check("S7e", "7", "peak-shape anomaly check",
          [(r"peak[_-]?shape|anomalous.*tail|asymmetry", "peaks/*.py")]),
    Check("S7f", "7", "reject if characteristic line absent",
          [(r"reject.*char|nuclide.*rejected|absent.*char", "identification/*.py")]),
    Check("S7g", "7", "additional isolated line + intensity-ratio test",
          [(r"intensity[_-]?ratio|S_i.*S_j|I_gamma_ratio", "identification/*.py")]),
    Check("S7h", "7", "chi2_intensity < 3 lenient / 1.5 strict",
          [(r"chi2.*<\s*3|chi2.*<\s*1\.5|intensity.*chi2", "identification/*.py")]),
    Check("S7i", "7", "library-directed search for I_gamma>=0.1%",
          [(r"0\.1\s*%|I_gamma.*0\.001|library[_-]?directed", "identification/*.py")]),
    Check("S7j", "7", "DPR equilibrium chains",
          [(r"DPR|chain[_-]?equilibrium|secular|equilibrium", "identification/*.py")]),
    Check("S7k", "7", "CI = log10(1/prod delta) per nuclide",
          [(r"confidence[_-]?index|\bCI\b|log10.*delta|log10.*prod", "identification/confidence.py")]),
    Check("S7l", "7", "DC = D_unident/D_total*100%",
          [(r"\bDC\b|dose[_-]?contribution|D_unident|D_total", "identification/completeness.py")]),

    Check("S8a", "8", "deconv ONLY after step 7 confirms list",
          [(r"after[_-]?identification|confirmed.*nuclide.*first|step\s*7", "peaks/deconvolve.py")]),
    Check("S8b", "8", "components fixed a-priori from confirmed nuclides",
          [(r"a[_-]?priori|expected[_-]?components|known[_-]?components", "peaks/deconvolve.py")]),
    Check("S8c", "8", "FWHM fixed from calibration",
          [(r"fwhm[_-]?fixed|fixed[_-]?fwhm|fwhm.*from.*cal", "peaks/deconvolve.py")]),
    Check("S8d", "8", "intensity ratios constrained to library; activity free",
          [(r"intensity[_-]?ratio.*constrain|library[_-]?ratio.*fix|activity.*free", "peaks/deconvolve.py")]),
    Check("S8e", "8", "HPGe Gauss+tail vs scint pure Gauss",
          [(r"gauss[_-]?plus[_-]?tail|low[_-]?energy[_-]?tail|hypermet|asymmetric.*gauss", "peaks/*.py")]),
    Check("S8f", "8", "STEP-AND-LINEAR continuum (Gilmore 9.7) — KEY GAP",
          [(r"step[_-]?linear|erfc|smooth[_-]?step|gilmore.*9\.7|step.*continuum", "peaks/*.py")]),
    Check("S8g", "8", "chi2/nu in [0.8, 1.5] QC gate",
          [(r"0\.8.*1\.5|chi2.*per.*dof|reduced.*chi2", "peaks/*.py")]),
    Check("S8h", "8", "area sigma from full covariance",
          [(r"covariance|pcov|full.*covariance|sigma.*cov", "peaks/*.py")]),
    Check("S8i", "8", "no free components beyond expected; chi2>1.5 -> step 7",
          [(r"return.*to.*step.*7|residual.*new.*candidate", "peaks/*.py")]),

    Check("S9a", "9", "per-FEP: ch E FWHM S cps sigma_cov significance",
          [(r"area.*sigma|cps.*sigma|S.*over.*sigma|significance", "peaks/area.py")]),
    Check("S9b", "9", "dead-time Lsrm Sec.15: t_m = A*Sum y + B*Sum y*i",
          [(r"dead[_-]?time.*A.*B|t_m\s*=|A\s*\*.*sum.*B\s*\*.*sum", "calibration/*.py")]),
    Check("S9c", "9", "flag activities uncorrected if A/B unavailable",
          [(r"dead[_-]?time.*uncorrected|A[_-]?B.*missing|flag.*uncorrected", "**/*.py")]),
    Check("S9d", "9", "MDA per ISO 11929 with L_C / L_D quadratic",
          [(r"ISO[_-]?11929|L_C|L_D|currie", "identification/mda.py")]),
    Check("S9e", "9", "decay-correction A0 = A*exp(ln2*dt/Thalf)",
          [(r"decay[_-]?correct|exp.*ln2|half[_-]?life.*correct|T_half", "**/*.py")]),

    Check("S10a", "10", "secondary-peak classifier",
          [(r"backscatter|single[_-]?escape|double[_-]?escape|compton[_-]?edge|pile[_-]?up|bremsstrahlung", "**/*.py")]),
    Check("S10b", "10", "intrinsic catalogue per detector",
          [(r"intrinsic|138La|K[_-]?escape|Ge[_-]?escape", "**/*.py")]),
    Check("S10c", "10", "BG carryover after independent BG calibration",
          [(r"background[_-]?carryover|bg[_-]?lines.*subtracted|recalibrated.*bg", "**/*.py")]),
    Check("S10d", "10", "XRF Z-by-Z catalog (Z 13-95)",
          [(r"xrf[_-]?catalog|K_alpha|K_beta|L_alpha|Z\s*=\s*\d", "data/xrf*.py")]),
    Check("S10e", "10", "critical overlap disambig (Am241/W, I125, U235/Th)",
          [(r"Am[-_]?241.*W|W.*K_alpha.*59|U[-_]?235.*Th|disambig", "identification/disambiguate.py")]),

    Check("S11a", "11", "report header",
          [(r"sample[_-]?id|live[_-]?time|dead[_-]?time|geometry|environment", "reporting/*.py")]),
    Check("S11b", "11", "detector justification block",
          [(r"detector[_-]?type|detector.*justific|resolution[_-]?check", "reporting/*.py")]),
    Check("S11c", "11", "calibration block: coeffs+sigma, residual plot, stored-vs-rebuilt",
          [(r"residual[_-]?plot|stored.*reused|rebuilt|cal[_-]?coeff", "reporting/*.py")]),
    Check("S11d", "11", "primary FEP table",
          [(r"primary[_-]?FEP|FEP[_-]?table|fep_table", "reporting/*.py")]),
    Check("S11e", "11", "secondary-peak table with type column",
          [(r"secondary[_-]?peak.*type|sec_peak.*table", "reporting/*.py")]),
    Check("S11f", "11", "elemental XRF table",
          [(r"xrf[_-]?table|fluorescence_shield|IC_parent|fluorescence_matrix", "reporting/*.py")]),
    Check("S11g", "11", "identified nuclides with CI, DPR, activity, decay-corrected",
          [(r"identified[_-]?nuclides.*CI|nuclide[_-]?table|activity_decay_corrected", "reporting/*.py")]),
    Check("S11h", "11", "unidentified + DC metric block",
          [(r"unident.*table|DC[_-]?metric|dose[_-]?contribution.*report", "reporting/*.py")]),
    Check("S11i", "11", "spectrum plot cps/log/no-smooth/labels",
          [(r"plot.*log|cps.*log|log.*cps|spectrum_plot", "reporting/*.py")]),
    Check("S11j", "11", "per-multiplet deconv plots",
          [(r"multiplet.*plot|deconv.*plot|fit[_-]?plot", "reporting/*.py")]),
    Check("S11k", "11", "MDA table for standard suite",
          [(r"mda[_-]?table|MDA.*suite|standard[_-]?nuclide.*MDA", "reporting/*.py")]),
    Check("S11l", "11", "diagnostic block",
          [(r"diagnostic[_-]?block|diagnostics.*report|shielding.*diagnostic", "reporting/*.py")]),
    Check("S11m", "11", "version history at end of report",
          [(r"version[_-]?history|changelog|report.*version", "reporting/*.py")]),
]


def search_canon(check):
    hits = []
    for pat, file_glob in check.patterns:
        try:
            rgx = re.compile(pat, re.IGNORECASE)
        except re.error as e:
            hits.append({"pattern": pat, "error": str(e), "files": []})
            continue
        glob_pat = file_glob if "/" in file_glob or file_glob.startswith("**") else f"**/{file_glob}"
        files_found = []
        for p in CANON.rglob("*.py"):
            rel = p.relative_to(CANON).as_posix()
            if not fnmatch(rel, glob_pat):
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="replace")
            except Exception:
                continue
            matches = []
            for i, line in enumerate(text.splitlines(), 1):
                if rgx.search(line):
                    matches.append({"line": i, "text": line.strip()[:160]})
                    if len(matches) >= 2:
                        break
            if matches:
                files_found.append({"file": rel, "matches": matches})
        hits.append({"pattern": pat, "file_glob": glob_pat, "files": files_found})
    return {"check_id": check.id, "step": check.step, "requirement": check.requirement, "patterns": hits}


def status_of(result):
    n_pat = len(result["patterns"])
    n_hit = sum(1 for p in result["patterns"] if p.get("files"))
    if n_hit == n_pat and n_pat > 0:
        return "PASS"
    if n_hit > 0:
        return "PARTIAL"
    return "MISS"


def main():
    print(f"# Skill vs canon audit -- {len(CHECKS)} checkpoints\n")
    print(f"Canonical tree: scripts/gamma/\n")

    out = {"checks": []}
    summary = {"PASS": 0, "PARTIAL": 0, "MISS": 0}
    by_step = {}
    for ch in CHECKS:
        r = search_canon(ch)
        r["status"] = status_of(r)
        out["checks"].append(r)
        summary[r["status"]] += 1
        by_step.setdefault(ch.step, []).append(r)

    step_order = ["principles", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]
    for step in step_order:
        if step not in by_step:
            continue
        print(f"## Step {step}\n")
        print("| ID | Status | Requirement | Hits |")
        print("|---|---|---|---|")
        for r in by_step[step]:
            hit_files = []
            for p in r["patterns"]:
                for f in p.get("files", []):
                    hit_files.append(f["file"])
            hit_str = ", ".join(sorted(set(hit_files))[:3]) if hit_files else "--"
            req = r["requirement"][:75]
            print(f"| {r['check_id']} | **{r['status']}** | {req} | `{hit_str}` |")
        print()

    print(f"## Summary\n")
    total = len(CHECKS)
    print(f"- **PASS**:    {summary['PASS']:>2} / {total}")
    print(f"- **PARTIAL**: {summary['PARTIAL']:>2} / {total}")
    print(f"- **MISS**:    {summary['MISS']:>2} / {total}")
    print()
    print("### MISS (no implementation found)\n")
    for r in out["checks"]:
        if r["status"] == "MISS":
            print(f"- **{r['check_id']}** (step {r['step']}): {r['requirement']}")
    print()
    print("### PARTIAL (only some patterns matched)\n")
    for r in out["checks"]:
        if r["status"] == "PARTIAL":
            missed = [p["pattern"] for p in r["patterns"] if not p.get("files")]
            print(f"- **{r['check_id']}** (step {r['step']}): {r['requirement']}")
            if missed:
                print(f"  - missing: `{missed[0]}`")

    json_path = ROOT / "scripts" / "audit" / "skill_vs_canon_audit_result.json"
    json_path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n_JSON dump_: scripts/audit/skill_vs_canon_audit_result.json")


if __name__ == "__main__":
    main()