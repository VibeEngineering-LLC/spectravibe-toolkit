"""
Convert LSRM NuclideMaster .lib to BecqMoni ROI XML format.
Input:  <LSRM_ROOT>\\Work\\BG\\Gamma-1S\\Data\\Библиотека радионуклидов\\Электроды.lib
Output: %APPDATA%\\BecqMoni\\config\\ROI\\Электроды.xml

BecquerelCoefficient=0 (no detector calibration).
ROI window: ±10% of peak energy.
Lines with intensity < 1% are skipped.
"""
import os
import sys
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

# Override via env: GAMMA_LIB_PATH / GAMMA_ROI_OUT
LIB_PATH = Path(os.environ.get(
    "GAMMA_LIB_PATH",
    r"C:\LSRM\Work\BG\Gamma-1S\Data\Библиотека радионуклидов\Электроды.lib"))
OUT_PATH = Path(os.environ.get(
    "GAMMA_ROI_OUT",
    os.path.join(os.environ.get("APPDATA", "."), r"BecqMoni\config\ROI\Электроды.xml")))

COLORS = {
    "W": "#FF8C00",
    "Th-232": "#0000FF",
}
DEFAULT_COLOR = "#FF0000"
MIN_INTENSITY = 1.0   # %
ROI_FRAC = 0.10       # ±10%


def half_life_years(nuc_elem):
    val = nuc_elem.get("half_life_value", "0").replace(",", ".")
    unit = nuc_elem.get("half_life_unit", "year")
    try:
        v = float(val)
    except ValueError:
        return 0.0
    factors = {"year": 1, "day": 1/365.25, "hour": 1/8766,
               "min": 1/525960, "sec": 1/3.156e7}
    return v * factors.get(unit, 1)


def roi_limits(e_kev, frac=ROI_FRAC):
    return round(e_kev * (1 - frac)), round(e_kev * (1 + frac))


def main():
    raw = LIB_PATH.read_bytes()
    text = raw.decode("windows-1251")
    # strip XML declaration so ElementTree doesn't choke on encoding attr
    if text.startswith("<?xml"):
        text = text[text.index("?>") + 2:].lstrip()
    root = ET.fromstring(text)

    now = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.0000000+03:00")
    rois = []

    for nuc in root.findall("Nuclide"):
        nuc_name = nuc.get("name")
        t_half = half_life_years(nuc)
        color = COLORS.get(nuc_name, DEFAULT_COLOR)

        for line in nuc.findall("Line"):
            energy = float(line.get("energy", "0").replace(",", "."))
            intensity = float(line.get("intensity", "0").replace(",", "."))
            line_type = line.get("line_type", "G")

            if intensity < MIN_INTENSITY:
                continue

            lo, hi = roi_limits(energy)

            if line_type == "X":
                xlt = line.get("xline_type", "")
                child = line.get("xray_child_nuclide", nuc_name)
                roi_name = f"{child} XRF {xlt} {energy:.2f}"
            else:
                dbid = line.get("dbid", "")
                parent = dbid.split("_")[0] if dbid else nuc_name
                roi_name = f"{parent} {energy:.3f} keV"

            note = (f"lib: {LIB_PATH.name} | nuclide: {nuc_name} | "
                    f"type: {line_type} | I={intensity:.2f}%")

            rois.append(dict(
                name=roi_name, enabled="true", peak=energy,
                lo=lo, hi=hi, color=color,
                bq_coef=0.0, bq_err=0.0,
                half_life=t_half, intensity=intensity,
                note=note,
            ))

    # Build ROI XML blocks
    blocks = []
    for r in rois:
        blocks.append(f"""\
    <ROIDefinitionData>
      <Name>{r["name"]}</Name>
      <Enabled>{r["enabled"]}</Enabled>
      <PeakEnergy>{r["peak"]:.3f}</PeakEnergy>
      <LowerLimit>{r["lo"]}</LowerLimit>
      <UpperLimit>{r["hi"]}</UpperLimit>
      <Color>{r["color"]}</Color>
      <BecquerelCoefficient>{r["bq_coef"]}</BecquerelCoefficient>
      <BecquerelCoefficientError>{r["bq_err"]}</BecquerelCoefficientError>
      <HalfLife>{r["half_life"]:.6g}</HalfLife>
      <Intencity>{r["intensity"]:.4f}</Intencity>
      <Note><![CDATA[{r["note"]}]]></Note>
      <ROIPrimitives>
        <ROISimpleDifferenceData>
          <PrimitiveType>BG difference</PrimitiveType>
          <OperationType>Addition</OperationType>
          <Coefficient>1</Coefficient>
          <CoefficientError>0</CoefficientError>
          <Note />
          <LowerLimit>{r["lo"]}</LowerLimit>
          <UpperLimit>{r["hi"]}</UpperLimit>
        </ROISimpleDifferenceData>
      </ROIPrimitives>
    </ROIDefinitionData>""")

    guid = str(uuid.uuid4())
    roi_block = "\n".join(blocks)

    xml_out = f"""<?xml version="1.0"?>
<ROIConfigData xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <FormatVersion>120920</FormatVersion>
  <Guid>{guid}</Guid>
  <Name>Электроды</Name>
  <LastUpdated>{now}</LastUpdated>
  <ROIDefinitions>
{roi_block}
  </ROIDefinitions>
  <ROIEfficiency />
  <Note />
</ROIConfigData>"""

    OUT_PATH.write_text(xml_out, encoding="utf-8")
    print(f"OK: {len(rois)} ROI -> {OUT_PATH}")
    for r in rois:
        print(f"  {r['name']:40s}  peak={r['peak']:8.3f}  lo={r['lo']:6d}  hi={r['hi']:6d}  I={r['intensity']:6.2f}%")


if __name__ == "__main__":
    main()
