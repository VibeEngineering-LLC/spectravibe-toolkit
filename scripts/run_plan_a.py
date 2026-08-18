"""Plan A wrapper (REPO-level) — Th-232 Marinelli canonical replay.

ЗАФИКСИРОВАНО ОПЕРАТОРОМ 2026-06-10: отчёты складывать в `demo_reports/` в отдельной
папке с датой и временем. Этот wrapper — точка входа Plan A с правильным дефолтным
output-dir. Под капотом — `scripts/run_skill.py` (F-398 / v1.18.28 production
orchestrator), который производит canonical-format bundle (включая главный
`sample_v2/<stem>_report.html` ~313 KB, Chart.js-based interactive viz).

ЗАФИКСИРОВАНО ОПЕРАТОРОМ 2026-06-10 (v1.2.7): после каждого Plan A прогона
автоматически вычисляется `cert_zcheck` — z-test + Δ% против сертификата для
каждого confirmed nuclide (Phase 1 exit-criteria: |Δ%| ≤ 15 AND |z| ≤ 2).
Defaults — для Th-232 Marinelli ОИСН-16 (1940 ± 6% Бк/кг).

Reference (форма отчёта):
    demo_reports/v1_18_32_th232_canonical/sample_v2/Th232_Маринелли_0cm_report.html

Usage (из корня репо):
    PYTHONIOENCODING=utf-8 python scripts/run_plan_a.py

Опциональные параметры через env vars:
    GAMMA_SAMPLE    — путь к sample .spe (default: Th-232 Marinelli ОИСН-16)
    GAMMA_BG        — путь к bg .spe
    GAMMA_MASS      — масса в кг (default 1.6)
    GAMMA_OUTPUT    — override output-dir (default: demo_reports/<ts>_<stem>/)
    GAMMA_LABEL     — suffix в имени output-папки (default: <sample stem>)
    GAMMA_CERT_A    — certificate activity, Бк/кг (default 1940 = Th-232 ОИСН-16)
    GAMMA_CERT_REL  — certificate σ relative (default 0.06 = ±6%)
    GAMMA_CERT_NAME — описание (default "Th-232 ОИСН-16 17-09-2007, m=1.6 кг")
    GAMMA_NO_ZCHECK — "1" чтобы выключить cert_zcheck (для проб без сертификата)
"""
from __future__ import annotations
import sys, os, subprocess, time, datetime, re
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

DEFAULT_SAMPLE = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Маринелли\Th232_420-7-17_Маринелли_0cm.spe"
DEFAULT_BG = r"C:\LSRM\Work\BG\Gamma-1S\Spe - поверки\Поверка 2024\Фон закр кр\Фон закр кр вода_13.spe"
DEFAULT_MASS_KG = 1.6
DEFAULT_CERT_A = 1940.0
DEFAULT_CERT_REL = 0.06
DEFAULT_CERT_NAME = "Th-232 ОИСН-16 17-09-2007, m=1.6 кг"

SAMPLE = os.environ.get("GAMMA_SAMPLE", DEFAULT_SAMPLE)
BG = os.environ.get("GAMMA_BG", DEFAULT_BG)
MASS_KG = float(os.environ.get("GAMMA_MASS", DEFAULT_MASS_KG))
CERT_A = float(os.environ.get("GAMMA_CERT_A", DEFAULT_CERT_A))
CERT_REL = float(os.environ.get("GAMMA_CERT_REL", DEFAULT_CERT_REL))
CERT_NAME = os.environ.get("GAMMA_CERT_NAME", DEFAULT_CERT_NAME)
SKIP_ZCHECK = os.environ.get("GAMMA_NO_ZCHECK", "") == "1"


def slug(name: str) -> str:
    s = re.sub(r"[^\w\-\.А-Яа-яЁё]+", "_", name)
    return s.strip("_")[:60]


def compute_output_dir() -> Path:
    override = os.environ.get("GAMMA_OUTPUT")
    if override:
        return Path(override).resolve()
    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    sample_stem = Path(SAMPLE).stem
    label = os.environ.get("GAMMA_LABEL") or sample_stem
    folder = f"{ts}_{slug(label)}"
    return (REPO / "demo_reports" / folder).resolve()


CHECK_ARTIFACTS_TEMPLATES = [
    ("sample_v2/{stem}_report.html", "V2 hybrid HTML — ГЛАВНЫЙ"),
    ("sample_v2/{stem}_report.json", "V2 hybrid JSON"),
    ("sample_v2/{stem}_report.md", "V2 hybrid Markdown"),
]

CLEANUP_PATHS_TEMPLATES = [
    "sample",
    "v2_compare",
    ".phases",
    "run_skill.log",
    "run_skill_summary.json",
    "index.html",
]


def detect_actual_stem(output: Path):
    candidates = sorted((output / "sample_v2").glob("*_report.html"))
    if not candidates:
        candidates = sorted((output / "sample").glob("*_report.html"))
    if not candidates:
        return None
    return candidates[0].name[:-len("_report.html")]


def log(msg):
    print(f"[plan-a] {msg}", flush=True)


def _on_rmtree_error(func, path, exc_info):
    import os, stat
    try:
        os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
        func(path)
    except Exception:
        pass


def cleanup_extras(output: Path) -> list:
    import shutil
    cleaned = []
    for rel in CLEANUP_PATHS_TEMPLATES:
        p = output / rel
        if not p.exists():
            continue
        try:
            if p.is_dir():
                shutil.rmtree(p, onerror=_on_rmtree_error)
            else:
                p.unlink()
            cleaned.append(rel)
        except Exception as e:
            log(f"WARN: cleanup {rel} failed: {e!r}")
    return cleaned


def run_cert_zcheck(output: Path) -> int:
    """Auto cert_zcheck. Returns helper exit code (0=ALL PASS, 3=any FAIL, 1/2=helper error)."""
    helper = REPO / "scripts" / "cert_zcheck.py"
    if not helper.exists():
        log(f"WARN: cert_zcheck.py not found at {helper}, skip")
        return -1
    cmd = [
        sys.executable, str(helper), str(output),
        "--cert-A", str(CERT_A),
        "--cert-rel", str(CERT_REL),
        "--cert-name", CERT_NAME,
    ]
    log(f"executing: scripts/cert_zcheck.py {output.name} (cert={CERT_NAME})")
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    r = subprocess.run(cmd, cwd=str(REPO), env=env, text=True, encoding="utf-8", errors="replace")
    return r.returncode



def inject_zcheck_into_html(output: Path, html_path: Path) -> bool:
    """Inject cert z-check table into HTML report before </body>."""
    import json as _j
    zp = output / "cert_zcheck.json"
    if not zp.exists() or not html_path.exists():
        return False
    try:
        d = _j.loads(zp.read_text(encoding="utf-8"))
        nu = d.get("nuclides", [])
        s = d.get("summary", {})
        cert = d.get("certificate", {})
        if not nu:
            return False

        def _row(r):
            v = r.get("verdict", "—")
            vc = "#2d7d2d" if v == "PASS" else ("#888" if "upper" in v.lower() else "#c0392b")
            A = f"{r['A_Bq_per_kg']:.1f}" if r.get("A_Bq_per_kg") is not None else "<DL"
            dp = f"{r['delta_pct']:+.2f}%" if r.get("delta_pct") is not None else "-"
            az = f"{r['abs_z']:.2f}" if r.get("abs_z") is not None else "-"
            return (f'<tr style="border-bottom:1px solid #eee">'
                    f'<td style="padding:4px 8px">{r["nuclide"]}</td>'
                    f'<td style="padding:4px 8px">{A}</td>'
                    f'<td style="padding:4px 8px">{dp}</td>'
                    f'<td style="padding:4px 8px">{az}</td>'
                    f'<td style="padding:4px 8px;color:{vc};font-weight:bold">{v}</td></tr>\n')

        ap = s.get("all_pass", False)
        scol = "#2d7d2d" if ap else "#c0392b"
        sl = f'{s.get("n_passed")}/{s.get("n_confirmed")} PASS'
        sl += " — ALL PASS" if ap else " — FAIL"
        cn = cert.get("name", "")
        rows = "".join(_row(r) for r in nu)
        block = (
            '\n<div id="zcheck-section" style="margin:24px auto;max-width:900px;'
            'padding:16px;border:1px solid #ddd;border-radius:6px;font-family:Arial,\'Segoe UI\',sans-serif;font-size:13px;">\n'
            '<h3 style="margin:0 0 8px">Z-test по сертификату (Phase 1)</h3>\n'
            f'<div style="color:#666;font-size:12px;margin-bottom:4px">{cn}</div>\n'
            f'<div style="color:{scol};font-weight:bold;margin-bottom:10px">{sl}</div>\n'
            '<table style="border-collapse:collapse;width:100%;font-size:12px;">\n'
            '<thead><tr style="background:#f5f5f5;text-align:left;">'
            '<th style="padding:4px 8px">Нуклид</th>'
            '<th style="padding:4px 8px">A, Бк/кг</th>'
            '<th style="padding:4px 8px">Δ%</th>'
            '<th style="padding:4px 8px">|z|</th>'
            '<th style="padding:4px 8px">Вердикт</th>'
            '</tr></thead>\n'
            f'<tbody>\n{rows}</tbody>\n'
            '</table>\n'
            '<div style="margin-top:8px;font-size:11px;color:#555;border-top:1px solid #e8e8e8;padding-top:6px;line-height:1.65">'
            '<b>Δ%</b> = (A<sub>изм</sub>&nbsp;&minus;&nbsp;A<sub>серт</sub>)/A<sub>серт</sub>&nbsp;&times;&nbsp;100% &mdash; '
            'относительное отклонение измеренной активности от сертифицированного значения.<br>'
            '<b>|z|</b> = (A<sub>изм</sub>&nbsp;&minus;&nbsp;A<sub>серт</sub>)/&radic;(σ²<sub>изм</sub>&nbsp;+&nbsp;σ²<sub>серт</sub>) &mdash; '
            'нормированное отклонение с учётом неопределённостей измерения и сертификата.<br>'
            'Пороги Phase&nbsp;1: <b>|Δ%|&nbsp;&le;&nbsp;15%</b> &mdash; допустимое отклонение;&nbsp;&nbsp;'
            '<b>|z|&nbsp;&le;&nbsp;2</b> &mdash; статистическая совместимость на уровне 2σ.'
            '</div>\n'
            '</div>\n'
        )
        html = html_path.read_text(encoding="utf-8")
        if "</body>" not in html:
            return False
        html_path.write_text(html.replace("</body>", block + "</body>", 1), encoding="utf-8")
        return True
    except Exception as e:
        log(f"WARN: inject_zcheck failed: {e!r}")
        return False


def main():
    if not Path(SAMPLE).exists():
        log(f"FAIL: sample not found: {SAMPLE}")
        return 2
    if not Path(BG).exists():
        log(f"FAIL: bg not found: {BG}")
        return 2

    output = compute_output_dir()
    output.mkdir(parents=True, exist_ok=True)
    sample_stem = Path(SAMPLE).stem

    log(f"sample:    {SAMPLE}")
    log(f"bg:        {BG}")
    log(f"mass:      {MASS_KG} kg")
    log(f"cert:      {CERT_NAME} (A={CERT_A} ±{CERT_REL*100:.1f}%)" + (" [SKIPPED]" if SKIP_ZCHECK else ""))
    log(f"output:    {output}")
    log(f"stem:      {sample_stem}")

    runskill = REPO / "scripts" / "run_skill.py"
    if not runskill.exists():
        log(f"FAIL: run_skill.py not found at {runskill}")
        return 2

    cmd = [
        sys.executable, str(runskill),
        SAMPLE,
        "--background", BG,
        "--mass", str(MASS_KG),
        "--include-v2",
        "--output-dir", str(output),
    ]
    log(f"executing: scripts/run_skill.py ... --output-dir {output.name}")
    t0 = time.time()
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    r = subprocess.run(cmd, cwd=str(REPO), env=env, text=True, encoding="utf-8", errors="replace")
    dt = time.time() - t0
    if r.returncode != 0:
        log(f"FAIL run_skill.py rc={r.returncode} in {dt:.1f}s")
        return 1
    log(f"run_skill.py OK in {dt:.1f}s")

    actual_stem = detect_actual_stem(output) or sample_stem
    if actual_stem != sample_stem:
        log(f"actual stem: {actual_stem} (filename stem was {sample_stem})")

    artifacts = []
    missing = []
    for rel_tpl, label in CHECK_ARTIFACTS_TEMPLATES:
        rel = rel_tpl.format(stem=actual_stem)
        p = output / rel
        if not p.exists():
            missing.append((rel, label))
            continue
        artifacts.append((rel, label, p))
        log(f"artifact OK: {rel} ({p.stat().st_size} bytes) — {label}")
    if missing:
        log("FAIL — отсутствуют артефакты:")
        for rel, label in missing:
            log(f"  - {rel} ({label})")
        return 2

    zcheck_rc = -1
    if not SKIP_ZCHECK:
        print()
        print("-" * 72)
        print("CERT Z-CHECK (Phase 1 exit-criteria: |Δ%| ≤ 15 AND |z| ≤ 2):")
        print("-" * 72)
        zcheck_rc = run_cert_zcheck(output)
        if zcheck_rc == 0:
            log("cert_zcheck: ALL PASS ✓")
        elif zcheck_rc == 3:
            log("cert_zcheck: FAIL (некоторые nuclides вне Phase 1 criteria)")
        elif zcheck_rc == -1:
            log("cert_zcheck: skipped (helper not found)")
        else:
            log(f"cert_zcheck: error rc={zcheck_rc}")
        main_html_path = output / "sample_v2" / f"{actual_stem}_report.html"
        if inject_zcheck_into_html(output, main_html_path):
            log("cert_zcheck: z-check table injected into HTML report")

    cleaned = cleanup_extras(output)
    if cleaned:
        log(f"cleanup: удалено {len(cleaned)} лишних путей (HARD-LOCK v1.2.8: отчёт = только sample_v2/)")

    print()
    print("=" * 72)
    print("АРТЕФАКТЫ ПРОГОНА (абсолютные пути для оператора):")
    print("=" * 72)
    for rel, label, p in artifacts:
        print(f"  • {label}:")
        print(f"      {p}")
    if not SKIP_ZCHECK:
        zpath = output / "cert_zcheck.json"
        if zpath.exists():
            print(f"  • cert_zcheck JSON:")
            print(f"      {zpath}")
    print("=" * 72)
    print("ГЛАВНЫЙ ОТЧЁТ (sample_v2 V2 hybrid HTML) — открыть в браузере:")
    main_html = output / "sample_v2" / f"{actual_stem}_report.html"
    print(f"  {main_html}")
    print("=" * 72)
    log("DONE")
    return 0 if (SKIP_ZCHECK or zcheck_rc in (0, -1)) else 4


if __name__ == "__main__":
    sys.exit(main())