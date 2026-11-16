"""Тесты профилей детекторов калибровки по форме собственных структур (AtomNano 5 PRO, AS2PRO).
Данные AS2PRO -- вне репозитория, только из переменных окружения; без них тесты AS2PRO пропускаются с причиной."""

import os
import glob
import json
import pathlib

import pytest

from gamma.calibration.intrinsic_shape import profiles
from gamma.calibration.intrinsic_shape import autopreset as ap
from gamma.calibration.intrinsic_shape import register as reg
from gamma.calibration.intrinsic_shape.core import load_spectrum, poly_eval
from gamma.calibration.intrinsic_shape.cli import main as cli_main
from gamma.io.readers import read_spectrum

# --- Константы ---

FIX = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "intrinsic_shape"
AN_BG = FIX / "labr3_background.xml"
AN_SAMPLE = FIX / "labr3_lu176_sample.xml"
AN8_BG = FIX / "labr3_8192_background.xml"
AN8_SAMPLE = FIX / "labr3_8192_sample_th_u.xml"

PRESETS = pathlib.Path(profiles.__file__).resolve().parent / "presets"

BASE_2048 = {"sample_a37": 28.657615750126208, "a37": 28.91276265432917,
             "sum": 846.9099397510021, "edge": 497.63476499112585}
BASE_8192 = {"sample_a37": 214.1551452912842, "a37": 117.69038532324872,
             "sum": 3404.3500554252264, "edge": 2001.3674735701677}

# 37 ± 6 кэВ (гамма-эквивалент 35.5–37.5 кэВ плюс ошибка калибровки файла ниже 100 кэВ), сумм-пик и β-край ± 3 %
E_TOL = {"a37": (37.0, 6.0), "sum": (1471.0, 44.0), "edge": (788.74, 24.0)}

# Данные AS2PRO (вне репозитория)
AS2PRO_DIR = os.environ.get("IS_AS2PRO_DIR", "")
AM_BG = os.environ.get("IS_AS2PRO_AM_BG", "")
AM_SAMPLE = os.environ.get("IS_AS2PRO_AM_SAMPLE", "")
# подмена профиля для мутационной проверки
AS2PRO_ID = os.environ.get("IS_TEST_AS2PRO_PROFILE", "labr_as2pro")


need_dir = pytest.mark.skipif(
    not (AS2PRO_DIR and os.path.isdir(AS2PRO_DIR)),
    reason="нет данных AS2PRO: задайте IS_AS2PRO_DIR (папка с фонами BG_*.xml; данные вне репозитория)",
)
need_am = pytest.mark.skipif(
    not (AM_BG and AM_SAMPLE and os.path.isfile(AM_BG) and os.path.isfile(AM_SAMPLE)),
    reason="нет пары AS2PRO с Am-241: задайте IS_AS2PRO_AM_BG и IS_AS2PRO_AM_SAMPLE (данные вне репозитория)",
)

# --- Вспомогательные функции ---

def _e_ok(cal, name, x):
    c, t = E_TOL[name]
    return abs(poly_eval(list(cal), x) - c) <= t


def _auto(bg, sample, pid):
    b = load_spectrum(bg)
    B = b["counts"]
    S = load_spectrum(sample, n_ch=len(B))["counts"]
    B[-1] = 0.0
    S[-1] = 0.0
    with profiles.using(pid):
        a = ap.autopreset(B, S, list(b["cal"]), verbose=False)
    return a, b["cal"]


# --- Фикстура ---

@pytest.fixture(autouse=True)
def _reset_profile():
    profiles.activate(None)
    yield
    profiles.activate(None)


# --- Тесты ---

def test_profiles_listed_and_default():
    assert set(profiles.ids()) >= {"labr_atomnano5pro", "labr_as2pro"}
    assert profiles.active_id() == "labr_atomnano5pro"
    for p in profiles.list_profiles():
        assert "method" in p and "presets" in p
        for kind in profiles.KINDS:
            assert profiles.preset_path(kind, p["id"]).exists()


def test_atomnano_profile_uses_original_presets():
    names = {"reference": "labr3_reference.json", "ecal": "labr3_ecal.json",
             "lowE": "labr3_lowE.json", "struct_windows": "labr3_struct_windows.json",
             "selftest": "labr3_selftest.json"}
    for k, fn in names.items():
        assert profiles.preset_path(k, "labr_atomnano5pro").name == fn
    assert ap.load_ref() == ap.load_ref(PRESETS / "labr3_reference.json")
    assert profiles.get("labr_atomnano5pro").get("sample_a37_near", False) is False
    assert profiles.get("labr_atomnano5pro").get("scan_tol_ch") == 3.0


def test_detect_by_device_name():
    assert profiles.detect("AtomSpectra Pro LaBr3(Ce)")[0] == "labr_as2pro"
    assert profiles.detect("AS2 PRO LaBr3")[0] == "labr_as2pro"   # второй шаблон match
    assert profiles.detect("AtomNano 5 PRO LaBr 2048")[0] == "labr_atomnano5pro"
    pid, msg = profiles.detect("Гамма-1С NaI 63x63")
    assert pid is None and "не распознан" in msg
    assert profiles.detect("")[0] is None


def test_resolve_unknown_device_warns():
    pid, w = profiles.resolve("auto", "неизвестный прибор")
    assert pid == profiles.DEFAULT_ID and len(w) == 1 and "по умолчанию" in w[0]
    pid, w = profiles.resolve("labr_as2pro", "AtomNano 5 PRO LaBr 2048")
    assert pid == "labr_as2pro"
    assert any("распознан как 'labr_atomnano5pro'" in s for s in w)


def test_resolve_unknown_profile_raises():
    with pytest.raises(KeyError):
        profiles.resolve("нет_такого")


def test_using_restores_previous():
    profiles.activate("labr_as2pro")
    with profiles.using("labr_atomnano5pro"):
        assert profiles.active_id() == "labr_atomnano5pro"
    assert profiles.active_id() == "labr_as2pro"


def test_require_method():
    profiles.require_method("la138_intrinsic", "labr_as2pro")
    with pytest.raises(ValueError):
        profiles.require_method("nai_external", "labr_as2pro")


def test_atomnano_fixtures_match_base():
    """Поведение AtomNano на обеих парах фикстур -- как в 2ca4d15 (значения базы)."""
    for bg, sm, base in ((AN_BG, AN_SAMPLE, BASE_2048), (AN8_BG, AN8_SAMPLE, BASE_8192)):
        dev = read_spectrum(str(bg)).detector_id
        pid, _ = profiles.resolve("auto", dev)
        assert pid == "labr_atomnano5pro"
        a, cal = _auto(bg, sm, pid)
        assert a["sample_a37"]["x"] == pytest.approx(base["sample_a37"], abs=1e-6)
        for nm in ("a37", "sum", "edge"):
            assert a["anchors"][nm]["x"] == pytest.approx(base[nm], abs=1e-6)


@need_dir
def test_as2pro_backgrounds_anchors():
    """Все фоны AS2PRO папки (не независимая выборка: допуски подбирались на них же): три опоры уверенные И их энергия
    по калибровке файла в допуске E_TOL; ловит допуск scan_tol_ch 3 вместо 12 и подмену профиля."""
    files = sorted(glob.glob(os.path.join(AS2PRO_DIR, "BG_*.xml")))
    assert files, AS2PRO_DIR
    for path in files:
        h = os.path.basename(path)
        a, cal = _auto(path, path, AS2PRO_ID)
        bad = [nm for nm in ("a37", "sum", "edge")
               if not (a["anchors"][nm]["ok"] and _e_ok(cal, nm, a["anchors"][nm]["x"]))]
        assert bad == [], (h, bad, a["warnings"])
        assert not any("расходится" in w for w in a["warnings"]), (h, a["warnings"])


@need_am
def test_as2pro_sample_a37_near_background_peak():
    """Образец с Am-241: пик ~37 кэВ образца берётся около пика фона, а не линия 59.5 кэВ (окно near профиля)."""
    a, cal = _auto(AM_BG, AM_SAMPLE, AS2PRO_ID)
    xb, xs = a["anchors"]["a37"]["x"], a["sample_a37"]["x"]
    assert abs(xs / xb - 1.0) < 0.15, (xb, xs)
    cs = load_spectrum(AM_SAMPLE)["cal"]
    c, t = E_TOL["a37"]
    assert abs(poly_eval(list(cs), xs) - c) <= t, poly_eval(list(cs), xs)


def test_legacy_profile_ids_are_aliases(tmp_path, capsys):
    """Устаревшие id (до переименования #IS-4) -- тот же профиль и предупреждение «устаревшее имя»; CLI их принимает."""
    for old, new in (("atomnano5pro_labr3", "labr_atomnano5pro"), ("as2pro_labr3", "labr_as2pro")):
        assert profiles.get(old)["id"] == new
        pid, w = profiles.resolve(old)
        assert pid == new and any("устаревшее имя" in s for s in w), w
        assert profiles.activate(old) == new and profiles.active_id() == new
    profiles.activate(None)
    rc = cli_main(["register", "--bg", str(AN_BG), "--sample", str(AN_SAMPLE), "--detector", "atomnano5pro_labr3", "--n-ms", "2", "--out", str(tmp_path / "r.json")])
    assert rc == 0 and profiles.active_id() == "labr_atomnano5pro"
    assert "устаревшее имя" in capsys.readouterr().out

def test_cli_detector_activates_profile(tmp_path):
    """CLI: --detector включает выбранный профиль, auto -- по имени прибора."""
    rc = cli_main(["register", "--bg", str(AN_BG), "--sample", str(AN_SAMPLE),
                   "--detector", "labr_as2pro", "--n-ms", "2", "--out", str(tmp_path / "r1.json")])
    assert rc == 0 and profiles.active_id() == "labr_as2pro"
    assert json.loads((tmp_path / "r1.json").read_text("utf-8"))["auto"]["identity"] is False

    profiles.activate(None)
    rc = cli_main(["register", "--bg", str(AN_BG), "--sample", str(AN_SAMPLE),
                   "--detector", "auto", "--n-ms", "2", "--out", str(tmp_path / "r2.json")])
    assert rc == 0 and profiles.active_id() == "labr_atomnano5pro"
    assert json.loads((tmp_path / "r2.json").read_text("utf-8"))["auto"]["identity"] is True


def test_cli_anchor_error_is_rc2(tmp_path, monkeypatch, capsys):
    """Пик-якорь не подгоняется -> сообщение и код 2, не трассировка."""
    def _boom(*a, **k):
        raise reg.AnchorFitError("якорь (пик ~37 кэВ) в образце не подгоняется (тест)")
    monkeypatch.setattr(reg, "anchor37", _boom)
    rc = cli_main(["register", "--bg", str(AN_BG), "--sample", str(AN_SAMPLE),
                   "--out", str(tmp_path / "r.json")])
    assert rc == 2
    assert "не подгоняется" in capsys.readouterr().err
