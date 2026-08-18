"""
Import X-ray lines from LSRM NuclideMaster X-ray.mdb into SQLite.
Source: C:\\LSRM\\NuclideMaster\\XRay_Viewer\\X-ray.mdb
Output: data/xray_lines_lsrm.sqlite

Run: python scripts/data/import_xray_mdb.py
"""
import csv
import sqlite3
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
CSV_PATH = PROJECT / "data" / "xray_lines_lsrm.csv"
DB_PATH = PROJECT / "data" / "xray_lines_lsrm.sqlite"

DDL = """
CREATE TABLE IF NOT EXISTS xray_lines (
    id        INTEGER PRIMARY KEY,
    element   TEXT    NOT NULL,
    Z         INTEGER NOT NULL,
    line      TEXT    NOT NULL,
    energy_keV REAL   NOT NULL,
    intensity REAL,
    error     REAL
);
CREATE INDEX IF NOT EXISTS idx_xray_z        ON xray_lines(Z);
CREATE INDEX IF NOT EXISTS idx_xray_element  ON xray_lines(element);
CREATE INDEX IF NOT EXISTS idx_xray_energy   ON xray_lines(energy_keV);
CREATE INDEX IF NOT EXISTS idx_xray_line     ON xray_lines(line);

CREATE TABLE IF NOT EXISTS xray_meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""

INSERT = """
INSERT INTO xray_lines (id, element, Z, line, energy_keV, intensity, error)
VALUES (?, ?, ?, ?, ?, ?, ?)
"""

META = [
    ("source_file", r"C:\LSRM\NuclideMaster\XRay_Viewer\X-ray.mdb"),
    ("import_date", "2026-07-04"),
    ("z_range", "5..104"),
    ("line_types", "Ka1,Ka2,Ka3,Kb1,Kb2,Kb3,Kb4,Kb5,KO,KP,La1,La2,Lb1,Lb2,Lb3,Lb4,Lb5,Lb6,Lg1,Lg2,Lg3,Lg6,Ll,Ln"),
    ("row_count", "1766"),
]


def main():
    if not CSV_PATH.exists():
        print(f"CSV not found: {CSV_PATH}", file=sys.stderr)
        sys.exit(1)

    if DB_PATH.exists():
        DB_PATH.unlink()
        print(f"Removed existing: {DB_PATH}")

    con = sqlite3.connect(DB_PATH)
    con.executescript(DDL)

    rows = []
    with open(CSV_PATH, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((
                int(row["number"]),
                row["element"],
                int(row["Z"]),
                row["line"],
                float(row["energy_keV"]),
                float(row["intensity"]) if row["intensity"] else None,
                float(row["error"]) if row["error"] else None,
            ))

    con.executemany(INSERT, rows)
    con.executemany("INSERT OR REPLACE INTO xray_meta VALUES (?,?)", META)
    con.commit()

    n = con.execute("SELECT COUNT(*) FROM xray_lines").fetchone()[0]
    zmin, zmax = con.execute("SELECT MIN(Z), MAX(Z) FROM xray_lines").fetchone()
    con.close()

    print(f"OK: {n} rows, Z={zmin}..{zmax}")
    print(f"DB: {DB_PATH}")


if __name__ == "__main__":
    main()
