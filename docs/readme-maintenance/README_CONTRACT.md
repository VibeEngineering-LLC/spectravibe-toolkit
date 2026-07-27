# Landing-page contract — `readme_en.md` / `readme_ru.md`

> Contract for the bilingual **landing / start pages** of SpectraVibe.
> Owner agent: **Agent C — Documentation** (`.claude/agents/agent-c-docs.md`).
> Operational checklist: [`UPDATE_CHECKLIST.md`](UPDATE_CHECKLIST.md).
> Automation skill: [`.claude/skills/readme-sync/SKILL.md`](../../.claude/skills/readme-sync/SKILL.md).

---

## 1. What the start pages are

The project has **two** start pages (entry points), latte-style:

- `readme_en.md` — English landing.
- `readme_ru.md` — Russian landing.

They are mutually linked via a language switcher (`English | Русский`) and live
at the **repository root**. Image assets live in `graphics_readme/`
(`logo_en.jpg`, `logo_ru.jpg`, `icon.jpg`, `screenshot.jpg`).

These are **operator-facing** marketing/overview pages. They are NOT the
changelog and NOT the methodology reference.

## 2. Hard boundaries (never cross)

- **`README.md` is pipeline-owned. NEVER edit it from the landing process.**
  It carries the version marker (`> **vX.Y.Z…`** — read by
  `scripts/version_check.py`) and the full cross-version changelog. It is
  changed only by the normal release process (Agent C version bump), never to
  keep the landings in sync.
- **Never edit Python** (`scripts/**`, any `*.py`, including `version_check.py`)
  as part of landing maintenance.
- The landings **link** to `README.md` as the "Changelog / version history"
  (the "What's new" section is a short pointer, not a copy).

## 3. RU ↔ EN parity rule

`readme_en.md` and `readme_ru.md` MUST stay **content-equivalent**:

- Same section set, same order, same facts, same code samples.
- Same badges, same asset layout (EN uses `logo_en.jpg`, RU uses `logo_ru.jpg`).
- When one side changes, the other side changes in the **same commit**.
- Slogans are fixed: EN = "When the spectrum becomes clear",
  RU = "Когда спектр становится понятным".

## 4. Facts only (anti-hallucination)

Every non-trivial claim about features / commands / formats MUST be backed by a
concrete file (e.g. `scripts/gamma/cli.py`, `FORMAT_REGISTRY.md`,
`docs/RUN_SKILL.md`, `detectors/`). If a fact cannot be verified in-repo, write
a `TODO`, do not invent. This mirrors Agent C's `anti-hallucination` policy.

## 5. Update triggers

Refresh **both** landings (per [`UPDATE_CHECKLIST.md`](UPDATE_CHECKLIST.md))
when any of these happen:

- A **new feature** lands that changes what an operator can do (new CLI flag,
  new format, new report artefact, new detector folder, etc.).
- A **development cycle is merged** to master (Agent A + B merge → Agent C wrap).
- A **`SKILL_VERSION` bump** changes the version shown in the badges.

## 6. Version-marker drift is acceptable for landing work

The repo is updated continuously, so the version badge in the landings, the
`SKILL_VERSION`, the `README.md` marker and the latest git tag may legitimately
be **out of sync** at any given moment. This is **expected** and is **not a
blocker** for writing or updating `readme_*`.

Consequence: `python scripts/version_check.py --allow-no-tag` may exit non-zero
purely because of a `git tag != SKILL_VERSION` drift. For landing maintenance
this is **fine** — landing files do not touch `README.md`, `json_report.py`, or
tags, so they cannot make this check worse. Only treat a *new* mismatch that the
landing edit itself introduced (it never should) as a real problem.

## 7. License

The project is **MIT** (`LICENSE` at repo root). The landings' License section
links to `LICENSE`; keep the badge and section in sync if the license changes.
