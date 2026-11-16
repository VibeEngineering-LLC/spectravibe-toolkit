"""Синтетическая пара NaI (1024 кан, ~3 кэВ/кан) для тестов окна: фон 20000 с и образец 1000 с с линией Cs-137 и дрейфом
усиления 1.5 %. Спектры -- модель tests/calibration/test_nai_shape_pipeline.synth, не измерения. Запуск из корня репозитория:
python apps/intrinsic-shape-gui/tests/make_nai_fixtures.py (перезаписывает tests/fixtures/nai_shape/*.xml; seed фиксирован)."""
import sys, pathlib
TK = pathlib.Path(__file__).resolve().parents[3]
sys.path[0:0] = [str(TK / "scripts"), str(TK / "tests" / "calibration")]
from numpy.random import default_rng
from test_nai_shape_pipeline import synth
from test_nai_shape_api import DEVICE_XML
from gamma.calibration.nai_shape import surrogate
OUT = TK / "tests" / "fixtures" / "nai_shape"
rng = default_rng(40)


def write(name, counts, live):
    pts = "\n".join("          <DataPoint>%d</DataPoint>" % c for c in counts)
    t = DEVICE_XML.replace("{NOTE}", "        <Note />\n").replace("{POINTS}", pts).replace("<Name>device</Name>", "<Name>NaI synthetic 1024</Name>")
    (OUT / name).write_text(t.replace(">1000<", ">%g<" % live), encoding="utf-8")


OUT.mkdir(parents=True, exist_ok=True)
write("nai_background.xml", rng.poisson(synth(scale=20.0)), 20000)
write("nai_sample_cs137.xml", surrogate.drift_events(rng.poisson(synth(scale=10.0, extra=((661.7, 3.0e6),))), [3.0, 1.015, 0.0], 1024, rng), 1000)
