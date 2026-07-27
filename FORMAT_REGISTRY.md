# Gamma Spectrum File Format Registry

Inventory of spectrum file formats — what our converter supports natively, what
we read but don't write yet, and the long-term roadmap based on the SpecUtils
reference (sandialabs/SpecUtils, the canonical C++ library that backs
InterSpec, BecqMoni and Cambio).

The converter goal is **bidirectional N-to-N conversion across all formats**.
The current release implements 4 native formats; the registry is designed so
new formats slot in by registering reader/writer hooks in
`gamma.io.format_registry`.

## Status legend

- ✅ Native read+write (lossless round-trip on the fields we model)
- 📖 Native read only
- 🚧 Planned (extension reserved in registry, no implementation yet)
- ❌ Out of scope (proprietary binary with no public spec)

## Currently supported

| Format | Ext(s) | R | W | Module | Notes |
|---|---|---|---|---|---|
| LSRM SpectraLine **binary** `.spe` | `.spe` | ✅ | ✅ | `gamma.io.lsrm_spe` | CP-1251 KEY=VALUE header + uint32-LE binary counts block. SpecUtils calls this `SpectraLine`. |
| LSRM SpectraLine **ASCII** `.spe` | `.spe` | ✅ | ✅ | `gamma.io.lsrm_spe_text` | `$DATE_MEA / $MEAS_TIM / $DATA / $ENER_FIT / $MCA_CAL / $ENDRECORD` text export by ЛСРМ SpectraLine. Disambiguated from binary variant by content sniff (`$` vs `SPECTR=`). Per ЛСРМ spec: `$MCA_CAL` first line is **N = number of coefficients** (not polynomial degree), NO unit suffix on coefficient line; lowercase `e` in exponential notation; optional `$SPEC_ID/$SPEC_REM/$ROI/$SHAPE_CAL`. SpecUtils calls this `LsrmSpe` and treats it as distinct from `SpeIaea`. |
| AtomSpectra / BecqMoni XML | `.xml` | ✅ | ✅ | `gamma.io.atomspectra_xml` + `gamma.io.becqmoni_xml` | `ResultDataFile` schema. Both AtomSpectra and BecqMoni emit this; SpecUtils calls it `RadiaCode` parser. |
| ANSI/IEEE N42.42-2012 | `.n42`, `.xml` | ✅ | ✅ | `gamma.io.n42_2012` | Namespace `http://physics.nist.gov/N42/2011/N42`. Used by InterSpec/SpecUtils. |

## Roadmap (planned but not yet implemented)

These are exposed by SpecUtils and may be added when fixtures and user
demand justify the work. The registry has placeholder entries that raise
`NotImplementedError` until a module lands.

### Read+Write candidates

| Format | Ext(s) | SpecUtils name | Priority |
|---|---|---|---|
| GADRAS PCF | `.pcf` | `Pcf` | Med — common in US DOE workflows |
| ORTEC CHN | `.chn` | `Chn` | Med — legacy ORTEC MCA, fixed-layout binary |
| Canberra CNF | `.cnf` | `Cnf` | Low — proprietary structure, reverse-engineered |
| SPC (ASCII / binary int / binary float) | `.spc` | `Spc`, `SpcBinaryInt`, `SpcBinaryFloat`, `SpcAscii` | Med — ORTEC GammaVision native |
| TKA (text) | `.tka` | `Tka` | Low — minimal counts-only |
| Exploranium GR-130/135 | `.bin` | `ExploraniumGr130v0`, `ExploraniumGr135v2` | Low — handheld backpack RIID |
| ANSI N42-2006 | `.n42` | `N42_2006` | Med — legacy variant for older instruments |
| **IAEA / ORTEC GammaVision SPE** | `.spe` | `SpeIaea` | Low — close cousin of LSRM ASCII SPE; reader could be added with minor tweaks (keV unit suffix, mandatory $SHAPE_CAL, integer $MEAS_TIM) |

### Read-only candidates

| Format | Ext(s) | SpecUtils name | Notes |
|---|---|---|---|
| Tracs MPS | `.mps` | `TracsMps` | |
| Aram TXT/XML hybrid | `.txt`, `.xml` | `Aram` | |
| SPM Daily | `.txt` | `SPMDailyFile` | |
| Amptek MCA | `.mca` | `AmptekMca` | |
| MicroRaider XML | `.xml` | `MicroRaider` | |
| ORTEC list mode | `.lis` | `OrtecListMode` | Event-by-event timestamps |
| MultiAct | `.mca` | `MultiAct` | Partial support upstream |
| PHD | `.phd` | `Phd` | CTBTO format |
| LabZY LZS | `.lzs` | `Lzs` | |
| Scan data XML | `.xml` | `ScanDataXml` | |
| Bridgeport JSON | `.json` | `Json` | MCA-3000 |
| CAEN Hexagon GXML | `.gxml` | `CaenHexagonGXml` | |
| URI (QR-code) | n/a | `Uri` | Encoded as URL |
| TXT/CSV/TSV | `.txt`, `.csv`, `.tsv` | `TxtOrCsv` | Highly variable — sniffer-driven |

### Write-only output formats (no read counterpart in SpecUtils)

| Format | Ext(s) | SpecUtils name | Notes |
|---|---|---|---|
| Interactive D3 HTML | `.html` | `HtmlD3` | Diagnostics-only |
| Inja template | (variable) | `Template` | User-defined output |

## Architecture

Each format is one module with two public functions:

```python
def read_<format>(path: str, **opts) -> Spectrum: ...
def write_<format>(spec: Spectrum, path: str, **opts) -> None: ...
```

`gamma.io.format_registry` maintains:

- `EXT_TO_FORMAT: dict[str, str]` — file extension → format id
- `READERS: dict[str, Callable[[str], Spectrum]]`
- `WRITERS: dict[str, Callable[[Spectrum, str], None]]`
- `SNIFFERS: list[Callable[[bytes], str | None]]` — for ambiguous extensions
  (notably `.spe` LSRM binary vs LSRM ASCII, and `.xml` AtomSpectra/BecqMoni vs N42)

`gamma.io.convert.convert_spectrum(in_path, out_path, *, in_format=None, out_format=None)`
auto-detects formats and runs `read → write`. The intermediate Spectrum
dataclass is the canonical neutral representation; what gets carried
between formats is exactly what the dataclass models. Format-specific
metadata that doesn't map to Spectrum fields lives in `spec.extras` and
is preserved opportunistically on output.

## Lossy-conversion caveats

Not every field round-trips through every format:

- LSRM `.spe` ↔ N42-42: LSRM's FWHM polynomial and PEAKS table have no N42
  counterpart; preserved in `extras` but not written to N42.
- AtomSpectra XML → LSRM ASCII SPE: SampleInfo, DeviceConfigReference and the
  full FWHM calibration block are lost.
- N42-42-2012 → N42-2006: cross-references and instrument metadata
  hierarchy collapse.
- Any format → LSRM ASCII SPE: at most one polynomial calibration (low-to-high)
  is preserved.

The converter prints a one-line summary of dropped fields when
`--verbose` is passed.

## References

- SpecUtils source: <https://github.com/sandialabs/SpecUtils>
- InterSpec: <https://github.com/sandialabs/InterSpec>
- BecqMoni (Nuclear edition): <https://github.com/Am6er/BecqMoni>
- ANSI/IEEE N42.42-2012 schema: <https://www.nist.gov/programs-projects/ansiieee-n4242-standard>
- LSRM SpectraLine ASCII SPE format: ЛСРМ SpectraLine User Manual (RU),
  section "Экспорт спектра в текстовом формате". Shares `$`-section
  layout with the ORTEC GammaVision / MAESTRO IAEA SPE format but
  uses LSRM-specific conventions ($MCA_CAL first line = N, no `keV`
  suffix, lowercase exponents).
