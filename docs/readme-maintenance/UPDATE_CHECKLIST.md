# Landing update checklist — `readme_en.md` / `readme_ru.md`

> Step-by-step procedure to keep both start pages current.
> Contract: [`README_CONTRACT.md`](README_CONTRACT.md). Run this whenever an
> update trigger fires (new feature / cycle merge / `SKILL_VERSION` bump).

---

## When to run

- A new operator-visible feature landed (CLI flag, format, report artefact,
  detector folder, …).
- A development cycle was merged to master.
- A `SKILL_VERSION` bump changed the version badge.

## Steps

1. **Identify the operator-visible delta.** What changed that an operator would
   see? Ignore pure internals (refactors, tests) — those do not belong on a
   landing page.

2. **Verify each new fact in-repo** (anti-hallucination). Cite a concrete file:
   - CLI flags / commands → `scripts/gamma/cli.py`, `docs/RUN_SKILL.md`
   - formats → `FORMAT_REGISTRY.md`
   - detectors → `detectors/`
   - data/config → `data/`
   If unverifiable → write `TODO`, do not invent.

3. **Edit `readme_en.md`** — update the affected section(s) only.

4. **Mirror into `readme_ru.md`** — same change, content-equivalent, in the
   **same commit** (RU↔EN parity, §3 of the contract).

5. **Badges** — if the version changed, bump the `version-vX.Y.Z` badge in
   *both* files. (Drift vs `SKILL_VERSION` / git tag is acceptable — see §6.)

6. **Assets** — if a screenshot/logo changed, update files under
   `graphics_readme/` and keep `logo_en.jpg` / `logo_ru.jpg` paired.

7. **Do NOT touch** `README.md` or any `*.py`. The "What's new" section stays a
   short pointer to `README.md`.

8. **Safety check** (informational):

   ```bash
   python scripts/version_check.py --allow-no-tag
   ```

   A non-zero exit caused by `git tag != SKILL_VERSION` drift is **expected and
   acceptable** for landing-only edits (contract §6). Confirm only that your
   landing edit did not change `README.md` / `json_report.py`:

   ```bash
   git status --short -- README.md '*.py' scripts/
   ```

   This should list nothing.

9. **Commit** (conventional, no push unless explicitly asked):

   ```
   docs(readme): sync landing pages for <feature/version>
   ```

## Quick sanity grep

```bash
# Both files should have the same section headers, in the same order:
grep -nE '^## ' readme_en.md
grep -nE '^## ' readme_ru.md
```
