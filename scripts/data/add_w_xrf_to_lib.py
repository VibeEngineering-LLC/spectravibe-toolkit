"""
Add W (tungsten) XRF nuclide to Electrody.lib LSRM library file.
Source: data/xray_lines_lsrm.sqlite (Z=74, from LSRM NuclideMaster X-ray.mdb).
Encoding: windows-1251 (as original file).

Run: python scripts/data/add_w_xrf_to_lib.py
"""
import re
import sqlite3
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
LIB_PATH = Path(r"C:\LSRM\Work\BG\Gamma-1S\Data") / "Библиотека радионуклидов" / "Электроды.lib"
DB_PATH = PROJECT / "data" / "xray_lines_lsrm.sqlite"
ENC = "windows-1251"

NUCLIDE_NAME = "W"


def fmt(v: float) -> str:
    """Float -> string with comma decimal separator (LSRM convention)."""
    s = f"{v:g}"
    return s.replace(".", ",")


def xline_type_lsrm(line_name: str) -> str:
    """Convert line name (Ka1, Kb2, La1...) to LSRM xline_type (KA1, KB2, LA1...)."""
    return line_name.upper()


def build_nuclide_xml(rows: list) -> str:
    # W: stable, atomic_mass=183.84
    lines_xml = []
    for line_name, energy_keV, intensity, error in rows:
        xlt = xline_type_lsrm(line_name)
        dbid = f"W_{line_name}"
        d_intensity = error if error is not None else 0.0
        lines_xml.append(
            f'    <Line energy="{fmt(energy_keV)}" d_energy="0" '
            f'intensity="{fmt(intensity)}" d_intensity="{fmt(d_intensity)}" '
            f'line_type="X" dbid="{dbid}" xline_type="{xlt}" '
            f'xray_child_nuclide="W" source="NuclideMaster"/>'
        )

    lines_block = "\n".join(lines_xml)
    return (
        f'  <Nuclide name="W" half_life_value="9E+99" half_life_unit="year" '
        f'gamma_constant="0" atomic_mass="183,84" d_atomic_mass="0,01" '
        f'nuc_num="9074" with_daughter="false" source="NuclideMaster">\n'
        f'{lines_block}\n'
        f'  </Nuclide>'
    )


def main():
    con = sqlite3.connect(DB_PATH)
    rows = con.execute(
        "SELECT line, energy_keV, intensity, error FROM xray_lines "
        "WHERE Z=74 ORDER BY energy_keV"
    ).fetchall()
    con.close()
    print(f"W lines from DB: {len(rows)}")

    raw = LIB_PATH.read_bytes()
    text = raw.decode(ENC)

    if 'name="W"' in text:
        print("W nuclide already present — aborting")
        return

    nuclide_xml = build_nuclide_xml(rows)

    # Insert before closing </Library>
    close_tag = "</Library>"
    if close_tag not in text:
        print(f"ERROR: '{close_tag}' not found in file")
        return

    new_text = text.replace(close_tag, f"{nuclide_xml}\n{close_tag}")

    LIB_PATH.write_bytes(new_text.encode(ENC))
    print(f"Written {len(rows)} W XRF lines to: {LIB_PATH}")


if __name__ == "__main__":
    main()
