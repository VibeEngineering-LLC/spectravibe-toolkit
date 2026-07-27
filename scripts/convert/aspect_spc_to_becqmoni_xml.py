"""#CONV-2: cs137_test1.spc (ASPECT) -> BecqMoni XML (round-trip via gamma.io)."""
from __future__ import annotations
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from gamma.io.aspect_spc import read_aspect_spc
from gamma.io.becqmoni_xml import write_becqmoni_xml

SRC = REPO / "tests/data/aspect_spc/cs137_test1.spc"
DST_DIR = REPO / "1_Version/v1.32.0/converted/2026-07-03_aspect_cs137_test1"
DST_DIR.mkdir(parents=True, exist_ok=True)
DST_XML = DST_DIR / "cs137_test1_becqmoni.xml"

spec = read_aspect_spc(str(SRC))
if not spec.sample_id:
    spec.sample_id = "cs137_test1"
if not spec.detector_id:
    spec.detector_id = "ASPECT_NaI_25mm2_handheld"
spec.comments = (
    "#CONV-2 ASPECT .spc -> BecqMoni XML. Source: tests/data/aspect_spc/cs137_test1.spc "
    "(SpecUtils issue #47 corpus, Cs-137 test measurement, 2016-01-10). Uncalibrated fixture "
    "(vendor gain=0, offset=0). Round-trip via gamma.io.becqmoni_xml.write_becqmoni_xml."
)
write_becqmoni_xml(spec, str(DST_XML), pretty=True)

# Verification via round-trip read
from gamma.io.atomspectra_xml import read_atomspectra_xml
back = read_atomspectra_xml(str(DST_XML))
assert back.n_channels == spec.n_channels, (back.n_channels, spec.n_channels)
assert int(back.counts.sum()) == int(spec.counts.sum())
assert abs(back.live_time - spec.live_time) < 1e-3
assert abs(back.real_time - spec.real_time) < 1e-3

size = DST_XML.stat().st_size
print("SRC   :", SRC, "->", SRC.stat().st_size, "B")
print("DST   :", DST_XML, "->", size, "B")
print("N     :", spec.n_channels, "==", back.n_channels)
print("sum   :", int(spec.counts.sum()), "==", int(back.counts.sum()))
print("live  :", f"{spec.live_time:.3f}", "==", f"{back.live_time:.3f}")
print("real  :", f"{spec.real_time:.3f}", "==", f"{back.real_time:.3f}")
print("start :", spec.start_datetime, "-> read back:", back.start_datetime)
print("OK: BecqMoni XML round-trip identity")
