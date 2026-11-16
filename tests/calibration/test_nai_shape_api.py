import sys, json, pathlib, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np, pytest
from numpy.random import default_rng
from scipy.special import ndtr

from gamma.calibration.nai_shape import api, pipeline, peakreg, strongcheck, surrogate, validate
from gamma.io.atomspectra_xml import read_atomspectra_xml
from test_nai_shape_pipeline import synth, spec, PROFILE, bg_pair
from test_nai_shape_scope import one_line_pair

ROOT = pathlib.Path(__file__).resolve().parents[2]
PROF_DIR = ROOT / "scripts" / "gamma" / "calibration" / "nai_shape" / "profiles"


def strong_pair():
    rng = default_rng(40)
    bg = spec(rng.poisson(synth(scale=20.0)), 20000.0, "bg")
    truth_s = synth(scale=10.0, extra=((661.7, 3.0e6),))
    sm = spec(surrogate.drift_events(rng.poisson(truth_s), [3.0, 1.015, 0.0], 1024, rng), 1000.0, "sm")
    return bg, sm


def inject_line(counts, E0, area):
    E = 3.0 * np.arange(len(counts))
    sig = 0.07 * np.sqrt(661.7 * E0) / 2.3548
    return np.asarray(counts, dtype=float) + area * (ndtr((E + 1.5 - E0) / sig) - ndtr((E - 1.5 - E0) / sig))


def csv_row(path):
    return open(path, encoding="utf-8").read().splitlines()[1].split(",")


def test_profiles_values():
    """Проверка значений профилей nai_gamma1s и nai_gs2020."""
    for fname in ("nai_gamma1s.json", "nai_gs2020.json"):
        p = json.load(open(PROF_DIR / fname, encoding="utf-8"))
        assert p["shape_min_bg_fraction"] == 0.5
        assert p["status"] == "experimental"
        assert p["min_windows"] == 3
        assert p["live_check_tol"] == 0.0
        assert p["live_check_nsigma"] == 3.0
        assert p["live_check_sigma_max"] == 0.005
        assert p["live_check_min_counts"] == 400
        assert p["id"] == fname.removesuffix(".json")


def test_profile_aliases_warn():
    """Алиасы профилей эмитируют DeprecationWarning."""
    with pytest.warns(DeprecationWarning):
        p = pipeline.load_profile(PROF_DIR / "gamma1s.json")
    assert p["id"] == "nai_gamma1s"
    with pytest.warns(DeprecationWarning):
        q = api.load_profile("gs2020")
    assert q["id"] == "nai_gs2020"


def test_choose_mode_threshold():
    """Пороговая логика выбора режима peaks/shape."""
    assert api.choose_mode(0.12, 0.5) == "peaks"
    assert api.choose_mode(0.45, 0.5) == "peaks"
    assert api.choose_mode(0.55, 0.5) == "shape"
    assert api.choose_mode(0.12, 0.5, force="shape") == "shape"
    with pytest.raises(ValueError):
        api.choose_mode(0.12, 0.5, force="x")


def test_choose_mode_boundary_inclusive():
    """Граничное значение f == thr включает режим shape."""
    assert api.choose_mode(0.5, 0.5) == "shape"


def test_scale_status_requires_all_conditions():
    """scale_status требует выполнения всех условий одновременно."""
    assert api.scale_status("shape", 0.6, 0.5, True, True) == (True, [])
    bad_cases = [
        ("peaks", 0.6, 0.5, True, True),
        ("shape", 0.12, 0.5, True, True),
        ("shape", 0.45, 0.5, True, True),
        ("shape", 0.6, 0.5, False, True),
        ("shape", 0.6, 0.5, None, True),
        ("shape", 0.6, 0.5, True, False),
    ]
    for args in bad_cases:
        v, reasons = api.scale_status(*args)
        assert v is False
        assert len(reasons) == 1

    v, reasons = api.scale_status("shape", 0.6, 0.5, True, True, n_windows=2)
    assert v is False
    assert len(reasons) == 1
    assert "нечем проверить" in reasons[0]

    assert api.scale_status("shape", 0.6, 0.5, True, True, n_windows=3)[0] is True

    v, reasons = api.scale_status("shape", 0.6, 0.5, True, True, bgpair_ok=False)
    assert v is False
    assert "фон-фон" in reasons[0]


def test_scale_status_boundary_inclusive():
    """Граничное f == thr при прочих условиях даёт verified=True."""
    assert api.scale_status("shape", 0.5, 0.5, True, True)[0] is True


def test_apply_status_red_verdict_gives_flag():
    """Красный вердикт присваивает флаг и обнуляет ok."""
    r = {"mode_auto": "shape", "verdict": {"ok": False, "n": 3}, "arrays": {"net": np.zeros(4)}}
    api.apply_status(r, 0.9, 0.5, True)
    assert r["scale_verified"] is False
    assert r["flag"] == "шкала не проверена"
    assert r["verdict"]["ok"] is None
    assert r["verdict"]["ok_raw"] is False
    assert r["verdict"]["counted"] is False
    assert r["status"] == "experimental"


def test_apply_status_green_verdict_verified():
    """Зелёный вердикт подтверждает проверку шкалы."""
    r = {"mode_auto": "shape", "verdict": {"ok": True, "n": 3}, "arrays": {"net": np.zeros(4)}}
    api.apply_status(r, 0.9, 0.5, True)
    assert r["scale_verified"] is True
    assert r["flag"] == ""
    assert r["verdict"]["ok"] is True
    assert r["verdict"]["counted"] is True
    assert r["status"] == "experimental"
    assert r["thresholds"] == {"thr": 0.5, "min_windows": 3}   # ревизия 7: печать убрана, пороги для пересчёта


def test_bg_fraction_ignores_counts_above_range():
    """bg_fraction игнорирует каналы за пределами энергетического окна."""
    B = np.full(1024, 10.0)
    B[1000] += 1.0e6
    S = np.full(1024, 10.0)
    frac = peakreg.bg_fraction(spec(B, 100.0, "bg"), spec(S, 100.0, "sm"), PROFILE)
    assert abs(frac - 1.0) < 1e-9


def test_live_decision_rules():
    """Правила live_decision с пустым и заполненным профилем."""
    P = {}
    assert api.live_decision(1.0, 0.002, 1e5, 1e5, True, P) == (True, "")
    assert api.live_decision(1.02, 0.003, 1e5, 1e5, True, P)[0] is False
    ok, reason = api.live_decision(1.0, 0.002, 100, 1e5, True, P)
    assert ok is False
    assert "мало отсчётов" in reason
    ok, reason = api.live_decision(1.0, 0.006, 1e5, 1e5, True, P)
    assert ok is False
    assert "точность" in reason
    ok, reason = api.live_decision(1.0, 0.002, 1e5, 1e5, False, P)
    assert ok is False
    assert "заголовок" in reason

    P2 = PROFILE
    assert api.live_decision(1.0, 0.002, 1e5, 1e5, True, P2)[0] is True
    assert api.live_decision(1.02, 0.003, 1e5, 1e5, True, P2)[0] is False
    assert api.live_decision(1.0, 0.002, 100, 1e5, True, P2)[0] is False
    assert api.live_decision(1.0, 0.006, 1e5, 1e5, True, P2)[0] is False
    assert api.live_decision(1.0, 0.002, 1e5, 1e5, False, P2)[0] is False


def test_header_plausible():
    """Проверка правдоподобности заголовка live/real."""
    assert api.header_plausible({"live": 1001.0, "real": 1000.0}, 0.002) is False
    assert api.header_plausible({"live": 970.0, "real": 1000.0}, 0.002) is False
    assert api.header_plausible({"live": 999.6, "real": 1000.0}, 0.002) is True


def test_live_check_detects_inflated_header():
    """live_check выявляет накрученный live при совпадении real."""
    rng = default_rng(7)
    bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
    sm = spec(rng.poisson(synth(scale=10.0, extra=((661.7, 2.0e5),))), 500.0, "sm")
    a = api.live_check(bg, sm, PROFILE, ["Cs-137"])
    b = api.live_check(bg, dict(sm, live=515.0), PROFILE, ["Cs-137"])
    c = api.live_check(bg, dict(sm, live=515.0, real=515.0), PROFILE, ["Cs-137"])
    assert a["ok"] is True
    assert b["ok"] is False
    assert b["header_ok"] is False
    assert c["ok"] is False
    assert c["header_ok"] is True
    assert abs(c["f"] - a["f"]) / a["f"] < 0.02


def test_live_check_default_dead_fraction():
    """Мёртвая доля 3 % при пороге по умолчанию 0.002 — заголовок неправдоподобен (ручное добавление, мутант X11)."""
    rng = default_rng(7)
    bg = spec(rng.poisson(synth(scale=20.0)), 1000.0, "bg")
    sm = spec(rng.poisson(synth(scale=10.0)), 500.0, "sm")
    assert api.live_check(bg, dict(sm, real=500.0 / 0.97), PROFILE, ["K-40"])["header_ok"] is False


def test_run_background_kind_verified_for_true_background():
    """run(kind='background') на истинном фоне даёт verified."""
    bg, sm = bg_pair([2.0, 1.01, 0.0], 3)
    r = api.run(bg, sm, PROFILE, kind="background")
    assert r["kind"] == "background"
    assert r["verdict"]["rule"].startswith("z_dip")
    assert r["bgpair"]["ok"] is True
    assert r["verdict"]["ok_raw"] is True
    assert r["scale_verified"] is True
    assert r["verdict"]["counted"] is True
    assert "verdict_fit" in r


def test_run_background_with_line_is_flagged():
    """run(kind='background') с линией помечает как непроверенный."""
    bg, sm = bg_pair([2.0, 1.01, 0.0], 3)
    sm2 = dict(sm, counts=inject_line(sm["counts"], 661.7, 10000.0))
    r = api.run(bg, sm2, PROFILE, kind="background")
    assert r["bgpair"]["ok"] is False
    assert r["scale_verified"] is False
    assert r["flag"] == "шкала не проверена"
    assert any("фон-фон" in s for s in r["scale_reasons"])
    assert abs(r["bgpair"]["E_worst"] - 661.7) < 60


def test_background_verdict_needs_min_windows():
    """background_verdict требует минимальное число окон."""
    bg, sm = bg_pair([2.0, 1.01, 0.0], 3)
    r = api.run(bg, sm, PROFILE, kind="background")
    v0 = api.background_verdict(r, dict(PROFILE, null_windows=[]))
    assert v0["n"] == 0
    assert v0["ok"] is False
    v2 = api.background_verdict(r, dict(PROFILE, null_windows=PROFILE["null_windows"][:2]))
    assert v2["n"] == 2
    assert v2["ok"] is False


def test_run_rejects_bad_arguments():
    """run отклоняет некорректные аргументы."""
    bg, sm = bg_pair([0.0, 1.0, 0.0], 1)
    with pytest.raises(ValueError):
        api.run(bg, sm, PROFILE, kind="x")
    with pytest.raises(ValueError):
        api.run(bg, sm, PROFILE, kind="background", force="shape")
    with pytest.raises(ValueError):
        api.run(bg, sm, PROFILE, [], kind="sample")


def test_write_outputs_flag_and_status_in_all_files(tmp_path):
    """write_outputs записывает флаг и статус во все файлы."""
    bg, sm = strong_pair()
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    assert r["scale_verified"] is False
    p = api.write_outputs(r, sm, str(tmp_path), "t")

    csv_text = open(p["csv"], encoding="utf-8").read()
    lines = csv_text.splitlines()
    assert lines[0].endswith(",status,scale_verified,flag,verdict_ok,verdict_ok_raw")
    assert all(line.strip() for line in lines)
    row = csv_row(p["csv"])
    assert row[-5] == "experimental"
    assert row[-4] == "0"
    assert row[-3] == "шкала не проверена"

    j = json.load(open(p["json"], encoding="utf-8"))
    assert j["status"] == "experimental"
    assert j["flag"] == "шкала не проверена"
    assert j["scale_verified"] is False
    assert j["verdict"]["counted"] is False
    assert row[-2] == str(j["verdict"]["ok"])

    for xkey in ("xml_zeroed", "xml_signed"):
        s = read_atomspectra_xml(p[xkey])
        assert "[экспериментально]" in s.sample_id
        assert "[шкала не проверена]" in s.sample_id
        assert "экспериментально, шкала образца не проверяется" in s.comments
        assert len(s.counts) == 1024


def test_write_outputs_verified_still_experimental(tmp_path):
    """Даже verified-результат сохраняет статус experimental."""
    bg, sm = bg_pair([2.0, 1.01, 0.0], 3)
    r = api.run(bg, sm, PROFILE, kind="background")
    p = api.write_outputs(r, sm, str(tmp_path), "v")
    s = read_atomspectra_xml(p["xml_signed"])
    assert "[экспериментально]" in s.sample_id
    assert "[шкала не проверена]" not in s.sample_id
    assert "шкала проверена" in s.comments
    assert "экспериментально" in s.comments
    row = csv_row(p["csv"])
    assert row[-5] == "experimental"
    assert row[-4] == "1"
    assert row[-3] == ""
    assert json.load(open(p["json"], encoding="utf-8"))["status"] == "experimental"


def test_write_outputs_rejects_forged_result(tmp_path):
    """write_outputs отклоняет подделанные результаты."""
    bg, sm = strong_pair()
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    f1 = dict(r, scale_verified=True, flag="", scale_reasons=[], verdict=dict(r["verdict"], ok=True, counted=True))
    with pytest.raises(ValueError):
        api.write_outputs(f1, sm, str(tmp_path), "f1")
    f2 = dict(r, scale_reasons=r["scale_reasons"][1:])   # ревизия 7: причины не совпадают с пересчётом
    with pytest.raises(ValueError):
        api.write_outputs(f2, sm, str(tmp_path), "f2")
    f3 = {k: v for k, v in r.items() if k != "thresholds"}
    with pytest.raises(ValueError):
        api.write_outputs(f3, sm, str(tmp_path), "f3")


def test_internal_writers_refuse_or_flag(tmp_path):
    """Внутренние писатели ставят флаг или бросают RuntimeError."""
    bg, sm = bg_pair([0.0, 1.0, 0.0], 1)
    raw = pipeline._run_pair(bg, sm, PROFILE, mode="background", verbose=False)
    assert raw["scale_verified"] is False
    assert raw["flag"] == "шкала не проверена"
    assert raw["status"] == "experimental"
    with pytest.raises(ValueError):   # ревизия 7: без статуса api.run — отказ проверки согласованности
        pipeline._write_outputs(raw, sm, str(tmp_path), "x")
    nb = validate.null_pair(bg, sm, PROFILE)
    assert nb["scale_verified"] is False
    assert nb["flag"] == "шкала не проверена"


def test_register_by_lines_raw_is_flagged():
    """_register_by_lines возвращает непроверенный результат."""
    bg, sm = one_line_pair(5, [3.0, 1.0, 0.0])
    r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
    assert r["scale_verified"] is False
    assert r["flag"] == "шкала не проверена"
    assert r["status"] == "experimental"


# минимальная копия структуры файла GS2020: имя образца в CDATA
DEVICE_XML = """<?xml version="1.0" encoding="utf-8"?>
<ResultDataFile xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <FormatVersion>120920</FormatVersion>
  <ResultDataList>
    <ResultData>
      <SampleInfo>
        <Name><![CDATA[00000001]]></Name>
        <Location />
{NOTE}      </SampleInfo>
      <DeviceConfigReference>
        <Name>device</Name>
      </DeviceConfigReference>
      <StartTime>2026-01-01T00:00:00</StartTime>
      <EndTime>2026-01-01T00:16:40</EndTime>
      <PresetTime>1000</PresetTime>
      <EnergySpectrum>
        <NumberOfChannels>1024</NumberOfChannels>
        <ChannelPitch>1</ChannelPitch>
        <EnergyCalibration>
          <PolynomialOrder>2</PolynomialOrder>
          <Coefficients>
            <Coefficient>0</Coefficient>
            <Coefficient>3</Coefficient>
            <Coefficient>0</Coefficient>
          </Coefficients>
        </EnergyCalibration>
        <ValidPulseCount>0</ValidPulseCount>
        <TotalPulseCount>0</TotalPulseCount>
        <MeasurementTime>1000</MeasurementTime>
        <LiveTime>1000</LiveTime>
        <NumberOfSamples>1</NumberOfSamples>
        <Spectrum>
{POINTS}
        </Spectrum>
      </EnergySpectrum>
    </ResultData>
  </ResultDataList>
</ResultDataFile>
"""


def device_xml(tmp_path, with_note):
    note = "        <Note />\n" if with_note else ""
    points = "\n".join("          <DataPoint>0</DataPoint>" for _ in range(1024))
    p = tmp_path / ("device_note.xml" if with_note else "device_nonote.xml")
    p.write_text(DEVICE_XML.replace("{NOTE}", note).replace("{POINTS}", points), encoding="utf-8")
    return p


@pytest.mark.parametrize("with_note", [True, False])
def test_write_outputs_marks_device_xml_cdata(tmp_path, with_note):
    """Пометки CDATA-имени и наличия Note в выгруженном XML."""
    import xml.etree.ElementTree as ET
    bg, sm = strong_pair()
    src = device_xml(tmp_path, with_note)
    sm = dict(sm, path=str(src), n_read=1024)
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    out = tmp_path / "out"
    p = api.write_outputs(r, sm, str(out), "dev")
    for key in ("xml_zeroed", "xml_signed"):
        s = read_atomspectra_xml(p[key])
        assert "[экспериментально]" in s.sample_id
        assert "[шкала не проверена]" in s.sample_id
        assert s.sample_id.startswith("dev")
        assert "экспериментально, шкала образца не проверяется" in s.comments
        si = ET.parse(p[key]).getroot().find("ResultDataList/ResultData/SampleInfo")
        assert len(si.findall("Name")) == 1
        assert len(si.findall("Note")) == 1
        assert ET.parse(p[key]).getroot().find("ResultDataList/ResultData/StartTime") is not None


def test_write_outputs_refuses_when_mark_missing(tmp_path, monkeypatch):
    """Отказ на выгрузку, если _mark_sampleinfo не применяется."""
    bg, sm = strong_pair()
    src = device_xml(tmp_path, True)
    sm = dict(sm, path=str(src), n_read=1024)
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    monkeypatch.setattr(api, "_mark_sampleinfo", lambda text, name, note_add: text)
    with pytest.raises(RuntimeError):
        api.write_outputs(r, sm, str(tmp_path / "o"), "dev")


def test_write_outputs_accepts_saved_result(tmp_path):
    """Принятие результата, восстановленного из сохранения."""
    bg, sm = strong_pair()
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    r2 = dict(r)
    p = api.write_outputs(r2, sm, str(tmp_path), "s")
    assert json.load(open(p["json"], encoding="utf-8"))["flag"] == "шкала не проверена"


def test_write_outputs_rejects_inconsistent_ok(tmp_path):
    """Отказ при несовместимом verdict.ok или status."""
    bg, sm = strong_pair()
    r = api.run(bg, sm, PROFILE, ["Cs-137"])
    bad = dict(r, verdict=dict(r["verdict"], ok=True))
    with pytest.raises(ValueError):
        api.write_outputs(bad, sm, str(tmp_path), "b")
    bad2 = dict(r, status="release")
    with pytest.raises(ValueError):
        api.write_outputs(bad2, sm, str(tmp_path), "b2")
def test_strongcheck_never_counted():
    """strongcheck.evaluate_pair никогда не помечает counted=True."""
    bg, sm = one_line_pair(5, [3.0, 1.0, 0.0])
    r = peakreg._register_by_lines(bg, sm, PROFILE, ["Cs-137"])
    e = strongcheck.evaluate_pair(r, bg, sm, PROFILE, ["Cs-137"])
    assert e["counted"] is False
    assert e["status"] in ("прошла", "не прошла")
    e2 = strongcheck.evaluate_pair(r, bg, sm, PROFILE, ["Cs-137"], min_windows=99)
    assert e2["counted"] is False
    assert e2["status"] == "нечем проверить"


def test_nested_status():
    """validate.nested_status проверяет временное и поканальное вложение."""
    B = {"start": "2026-10-13 23:00:00", "real": 1000.0, "counts": np.full(8, 10.0)}
    A = {"start": "2026-10-13 23:00:00", "real": 400.0, "counts": np.full(8, 4.0)}
    assert validate.nested_status(B, A)[0] is True

    ok, msg = validate.nested_status(B, dict(A, start="2026-10-13 22:59:00"))
    assert ok is False
    assert msg.startswith("невложенные")

    ok, msg = validate.nested_status(B, dict(A, start="2026-10-13 23:10:00", real=900.0))
    assert ok is False
    assert msg.startswith("невложенные")

    ok, msg = validate.nested_status(B, dict(A, counts=np.array([4.0] * 7 + [11.0])))
    assert ok is False
    assert msg.startswith("невложенные")
