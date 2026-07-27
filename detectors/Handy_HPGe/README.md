# Handy_HPGe Detector Class

## 1. Taxonomy lock — Crystal vs Station (2026-06-05)

| Crystal-class     | Station-instance   | Resolution |
|-------------------|--------------------|------------|
| Handy_HPGe        | Point-15cm         | TBD        |

State: "Crystal-class source pin: audit/_rag/visual_templates/SCHEMA.md:10"

Canonical name table entry for Handy_HPGe.

## 2. Detector description

This is a handheld coaxial HPGe spectrometer.

«UNVERIFIED — requires confirmation against detector passport»

Crystal material: Coaxial HPGe

«UNVERIFIED — requires confirmation against detector passport»

Crystal dimensions: TBD — requires operator verification against detector passport.

## 3. Vessel canonicalization

Geometry tags applicable to this class:
- Point-15cm (39)

Visual-template shape is determined by the tuple (crystal_class × vessel_class × useful_sample_volume_ml). Constructive part (vessel + d_eff) is portable between stations of the same crystal-class. Operational part (useful_volume + placement_distance) is per-station.

## 4. Runtime conflict-resolution rule (user lock 2026-06-05)

Geometry parameters at measurement time may diverge from any nominal spec. The folder name / .spe header / passport / operator metadata are four independent provenance layers, and in general they can diverge.

| Provenance Layer            | Rank |
|-----------------------------|------|
| operator_explicit_metadata  | 1    |
| passport_pdf_block          | 2    |
| spe_header                  | 3    |
| folder_name_hint            | 4    |
| lsrm_source_default         | 5    |

Conflict handling: offline pipeline logs to __geometry_conflicts + sets __needs_operator_review. Online analyzer MUST prompt operator — no silent resolution.

## 5. Layout

```
detectors/Handy_HPGe/
├── README.md
├── certificates/           (.gitkeep — detector passports, source certificates)
├── efficiency/             (.gitkeep — efficiency curves per geometry)
├── references/             (.gitkeep — reference documents)
├── reference_spectra/      (.gitkeep — reference .spe files)
└── data/
    └── SPECTRA_MANIFEST.json   (metadata-only, populated by F-300/W3)
```

Note: NO raw_lsrm/ subdirectory — see §7.

## 7. Local-only raw exclusion (F-115 / F-150)

raw_lsrm/ is NOT created for this class.
All spectra are metadata-only in data/SPECTRA_MANIFEST.json.
F-115: no operator absolute paths are committed.
F-150: no binary .spe / .zip files are committed.
LSRM path hint for reference only: Work\Handy\Handy(HPGe)\Spe\  (operator-side only, F-115)

## 8. Isolation policy (F-83)

This folder is isolated from other detectors. Algorithms in scripts/gamma/ are shared, but data assets here are valid only for the Handy_HPGe class.

## 11. Cross-refs

F-83, F-115, F-150, F-153, F-155, F-256, F-265, detector taxonomy lock 2026-06-05.

Note: SPECTRA_INDEX records for this class: 39 (source: audit/_plans/F-300_W0_recon_report.md §2).

## 12. Last release

v1.26.0 — F-300 W2 initial skeleton creation.