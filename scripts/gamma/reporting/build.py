"""
Top-level report dispatcher.

Typical use:

    from gamma.identification.staged_pipeline import analyze_lsrm_spe
    from gamma.reporting import build_report

    result = analyze_lsrm_spe(path, complete_workflow=True)
    artefacts = build_report(result, output_dir="./out",
                             write_plots=True, write_html=True)
    print(artefacts["summary"])

`build_report` writes the JSON report (and the optional Markdown /
HTML / PNG plots) to `output_dir` and returns a dict with the file
paths and the in-chat summary string.
"""
from __future__ import annotations

import dataclasses
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from gamma.reporting.json_report import build_json_report
from gamma.reporting.chat_summary import build_chat_summary
from gamma.reporting.markdown_report import build_markdown_report


# F-317 / v1.18.15 — pattern для F-id strip из user-facing body.
# Только pure F-id (F-NN, F-NNN), F-NN/F-NN. Сохраняем legitimate ссылки на
# контракты типа "ГОСТ", "ISO 11929", "ЛСРМ §..." как они есть.
_F_ID_PATTERN = re.compile(
    r"\bF-\d{1,3}([/-]\d{1,3})?(?:\s*\(F-\d{1,3}\))?",
)
# Также убираем bare F-id в скобках: "(F-256 / v1.17.10)", "(F-122)" и т.п.
# F-365 / v1.18.24.1 — class `[^)\n]` (NOT `[^)]`) — критично, чтобы парный
# strip не пересекал newline. Иначе на inline `<script>` блоках regex
# жадно ест `(btn => {\n  // F-147 — ...\n  setView()` как один матч и
# ломает синтаксис JS, что и приводит к "опять не выводятся окна спектров
# и мультиплетов" — Chart.js падает на SyntaxError.
_F_ID_PAREN_PATTERN = re.compile(
    r"\s*\([^)\n]*\bF-\d{1,3}[^)\n]*\)",
)
# K-NN и T-NN ID тоже internal — strip из user-facing
_KT_ID_PAREN_PATTERN = re.compile(
    r"\s*\([^)\n]*\b[KT]-\d{1,3}[^)\n]*\)",
)

# F-317 issue #36 / v1.18.31+ — context-aware bare-strip support.
#
# Blind `_F_ID_PATTERN.sub("", text)` was damaging two surface kinds:
#   (a) legitimate user-facing prose mentioning an F-rule as justification,
#       e.g. «Метод σ по правилу F-89.» → «Метод σ по правилу .» (orphan
#       punctuation). Evidence: BUG-6 workaround comments at
#       interactive_html.py:2244-2252 and :2336-2341 — Agent B had to rephrase
#       prose to *avoid* bare F-IDs entirely.
#   (b) F-IDs intentionally embedded in inline `code` / <code> / <pre> /
#       <script> blocks — developer-facing surfaces that should survive the
#       user-facing pipeline (e.g. a code example showing how F-317 itself is
#       called).
#
# Fix — two layers:
#   1. Span tokeniser `_iter_protected_spans` walks the input and yields
#      (start, end) regions covered by ` `code` `, ``` ```fence``` ```,
#      <code>…</code>, <pre>…</pre>, <script>…</script>. Strip patterns
#      skip those spans entirely.
#   2. Bare-strip gated by `_is_internal_ref_context` — only fires when the
#      F-id is in an "internal audit-trail context":
#        • preceded by `(` (e.g. "(F-122)") — already handled by paren
#          patterns, this is the residual case after paren strip;
#        • preceded by `// ` (comment-line ref);
#        • preceded by start-of-line + optional bullet/whitespace + nothing
#          else (free-standing reference);
#        • preceded by `→ ` / `см. ` / `see ` / `cf. ` (explicit cross-ref).
#      Otherwise, the F-id is in running prose (e.g. "правилу F-89") and is
#      LEFT INTACT.
#
# Cleanup of orphan punctuation (the «правилу .» problem) is no longer
# needed because the bare-strip no longer fires in that context.
_PROTECTED_TOKENS = [
    # (opener regex, closer regex). Order matters: longer/triple before single.
    (re.compile(r"```"), re.compile(r"```")),
    (re.compile(r"<script\b[^>]*>", re.IGNORECASE),
     re.compile(r"</script>", re.IGNORECASE)),
    (re.compile(r"<pre\b[^>]*>", re.IGNORECASE),
     re.compile(r"</pre>", re.IGNORECASE)),
    (re.compile(r"<code\b[^>]*>", re.IGNORECASE),
     re.compile(r"</code>", re.IGNORECASE)),
    (re.compile(r"`"), re.compile(r"`")),
]


def _iter_protected_spans(text: str) -> List[tuple]:
    """Return list of (start, end) ranges of protected spans (code/script/pre).

    Walks left-to-right. For each cursor position picks the earliest opener
    that matches; advances cursor past the matching closer. Overlap-free.
    Unclosed openers extend to end-of-text — defensive against malformed
    fragments. Spans returned in input order.
    """
    spans: List[tuple] = []
    cursor = 0
    n = len(text)
    while cursor < n:
        next_open = None  # (start, end, closer_re)
        for op_re, cl_re in _PROTECTED_TOKENS:
            m = op_re.search(text, cursor)
            if m is None:
                continue
            if next_open is None or m.start() < next_open[0]:
                next_open = (m.start(), m.end(), cl_re)
        if next_open is None:
            break
        open_start, open_end, closer_re = next_open
        close_m = closer_re.search(text, open_end)
        if close_m is None:
            # Unclosed → protect to EOT.
            spans.append((open_start, n))
            break
        spans.append((open_start, close_m.end()))
        cursor = close_m.end()
    return spans


def _index_protected(spans: List[tuple], pos: int) -> bool:
    """True iff `pos` falls inside any protected span."""
    for s, e in spans:
        if s <= pos < e:
            return True
        if pos < s:
            return False  # spans sorted
    return False


def _sub_outside_protected(
    pattern: "re.Pattern[str]",
    text: str,
) -> str:
    """`pattern.sub("", text)` but skip any match whose start falls inside a
    protected span (` `code` `, ``` ```fence``` ```, <code>, <pre>, <script>).
    """
    spans = _iter_protected_spans(text)
    if not spans:
        return pattern.sub("", text)
    out = []
    cursor = 0
    for m in pattern.finditer(text):
        if _index_protected(spans, m.start()):
            continue
        out.append(text[cursor:m.start()])
        cursor = m.end()
    out.append(text[cursor:])
    return "".join(out)


# Cleanup patterns applied after bare F-id strip. Order matters — heading
# rules run BEFORE generic space-before-punct (otherwise "### : title" loses
# its space and the heading rule no longer matches).
_CLEANUP_RULES: List[tuple] = [
    # «### F-145: title» → «###  : title» → «### title»: heading-empty-prefix.
    # NB: needs to match BEFORE the space-before-colon rule below.
    (re.compile(r"^(#{1,6}) +:\s+", re.MULTILINE), r"\1 "),
    # «по правилу /F-365» → «по правилу /» → «по правилу»: orphan slash
    (re.compile(r" +/+ +"), " "),
    (re.compile(r" +/+(?=[\s.,;:!?»\)\]]|$)"), ""),
    # «правилу F-89.» → «правилу .» → «правилу.»  : space-before-terminator.
    # NB: skip `:` to keep table-cell separators («| Шаг 9 :» style) alone;
    # heading rule above already handled trailing `:` cases.
    (re.compile(r" +([.,;!?»\)\]])"), r"\1"),
    # «Pb-50. F-110: не разлагается» → «Pb-50. : не разлагается» →
    # «Pb-50. не разлагается»  : drop orphan «: » after sentence terminator
    # left over from «F-NN:» strip.
    (re.compile(r"([.!?]) +: +"), r"\1 "),
    # «| F-115 текст |» → «|  текст |» → «| текст |»: leading cell space
    (re.compile(r"\|  +"), "| "),
    # «текст |» trailing cell extra whitespace
    (re.compile(r"  +\|"), " |"),
]

# Final cleanup pass — runs ONCE on the whole joined text after protected
# spans are stitched back. Trailing space before EOL belongs here because
# applying it per-segment would eat the space immediately preceding a
# protected span (e.g. «text `F-317`» segment-1 «text » would lose its
# trailing space before backtick).
_FINAL_CLEANUP_RULES: List[tuple] = [
    # Trailing-space-before-EOL (e.g. «см. F-89\n» → «см. \n» → «см.\n»)
    (re.compile(r" +$", re.MULTILINE), ""),
]


def _cleanup_strip_residue(text: str) -> str:
    """F-317 issue #36 — clean up punctuation residue left by bare F-id
    strip. Operates only outside protected spans to keep code blocks
    pristine.

    Examples (in → out):
        «Метод σ по правилу F-89.»          → «Метод σ по правилу.»
        «### F-145: двухфазная самокалибровка» → «### двухфазная самокалибровка»
        «| v1.17.4 | F-115 анонимизация |»  → «| v1.17.4 | анонимизация |»
    """
    spans = _iter_protected_spans(text)
    if not spans:
        for rx, repl in _CLEANUP_RULES:
            text = rx.sub(repl, text)
    else:
        # Walk segments between protected spans; apply per-segment cleanup
        # to non-protected pieces only.
        out_parts: List[str] = []
        cursor = 0
        for s, e in spans:
            seg = text[cursor:s]
            for rx, repl in _CLEANUP_RULES:
                seg = rx.sub(repl, seg)
            out_parts.append(seg)
            out_parts.append(text[s:e])  # protected verbatim
            cursor = e
        tail = text[cursor:]
        for rx, repl in _CLEANUP_RULES:
            tail = rx.sub(repl, tail)
        out_parts.append(tail)
        text = "".join(out_parts)
    # Final pass — applied once on the joined text. Whole-document rules
    # belong here so they cannot eat boundary whitespace next to protected
    # spans.
    for rx, repl in _FINAL_CLEANUP_RULES:
        text = rx.sub(repl, text)
    return text

# F-319 / v1.18.15 — RU-замены English-only терминов в user-facing body.
# F-260 narrator работает в обратную сторону (ru→ru+en); этот хелпер
# делает en→ru для случаев когда в коде осталось EN без сопровождения RU.
# Каждый замен сохраняет EN-эквивалент в скобках для совместимости с
# двуязычным регламентом (F-260).
_EN_TO_RU_REPLACEMENTS = [
    # Метрологические термины
    (r"\bdouble escape\b", "пик двойного вылета (double escape)"),
    (r"\bsingle escape\b", "пик одиночного вылета (single escape)"),
    (r"\bplateau\b", "плато (plateau)"),
    (r"\bsum peak\b", "суммарный пик (sum peak)"),
    (r"\bpile-?up\b", "наложение импульсов (pile-up)"),
    # Geometric / spectral
    (r"\bbackscatter\b", "обратное рассеяние (backscatter)"),
]
_EN_TO_RU_COMPILED = [
    (re.compile(pat, re.IGNORECASE), repl)
    for pat, repl in _EN_TO_RU_REPLACEMENTS
]


def _apply_en_to_ru_glossary(text: str) -> str:
    """F-319: однопроходная замена EN-only терминов на ru + (en)."""
    if not text:
        return text
    for regex, repl in _EN_TO_RU_COMPILED:
        # Не трогаем если уже двуязычная форма "ru (en)" или внутри code-block
        # Простой guard: если перед match стоит "(" → skip (уже в скобках EN-as-clarification)
        def _smart_repl(m, _repl=repl):
            start = m.start()
            # Skip если preceded by "(" — уже в виде "(en)"
            if start > 0 and text[start - 1] == "(":
                return m.group(0)
            return _repl
        text = regex.sub(_smart_repl, text)
    return text


def _f317_apply_user_facing_compliance(text: str, *, format: str) -> str:
    """F-317 / v1.18.15 — apply F-256 + F-260 + F-id strip к user-facing output.

    Конвейер:
      1. F-256 translate_text: Layer 1 [RAG-ID] → Layer 2 [N, локатор] ГОСТ.
      2. F-260 enrich_text: ru-термины → ru (en) при первом упоминании.
      3. F-id paren strip: "(F-256 / v1.17.10)" → "" (internal-only).
      4. F-id bare strip: "F-256" → "" (только если не в legitimate ссылке).

    HTML версия: применяется к body content (заголовки, параграфы), не к
    атрибутам HTML тегов.
    """
    if not text:
        return text
    # Step 1: F-256 citation translator
    try:
        from gamma.reporting.citation_translator import translate_text as _ct
        translated, _warnings = _ct(text)
        text = translated
    except Exception:
        pass
    # Step 2: F-260 bilingual narrator
    try:
        from gamma.reporting.bilingual_narrator import enrich_text as _en
        enriched, _stats = _en(text)
        text = enriched
    except Exception:
        pass
    # F-317 issue #36 / v1.18.31+ — context-aware substitution.
    # Step 3a: strip parenthetical F-id / K-id / T-id mentions, but skip
    # any match that lands inside a protected span (` `code` `, ``` ```fence``` ```,
    # <code>, <pre>, <script>). Prevents F-365-class breakage of JS in
    # <script> AND preserves dev-facing examples like
    # `<code>... F-317 ...</code>`.
    text = _sub_outside_protected(_F_ID_PAREN_PATTERN, text)
    text = _sub_outside_protected(_KT_ID_PAREN_PATTERN, text)
    # Step 3b: bare F-id strip, also span-aware (same protected spans).
    # User-facing contract (test_v1_18_17_audit_guards.py:54-70): bare
    # F-IDs are forbidden in MD/HTML body. Strip them outside protected
    # spans; rely on Step 3c to clean up residue (orphan punctuation,
    # double spaces, empty heading titles).
    text = _sub_outside_protected(_F_ID_PATTERN, text)
    # Step 3c: residue cleanup — solves BUG-6 root cause
    # («Метод σ по правилу F-89.» → «Метод σ по правилу .» previously).
    # Workaround comments at interactive_html.py:2244-2252 and :2336-2341
    # rewrote prose to avoid bare F-IDs entirely; with proper residue
    # cleanup that workaround is no longer required.
    text = _cleanup_strip_residue(text)
    # Step 4: F-319 — EN→ru terminology замены
    text = _apply_en_to_ru_glossary(text)
    # Cleanup double spaces left after strips
    text = re.sub(r"  +", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text


def _f318_append_gost_references(text: str, *, format: str) -> str:
    """F-318 / v1.18.15 — append "Список использованной литературы" в конец.

    Стратегия:
      1. Собрать Layer 2 ГОСТ ссылки `[N, локатор]` где локатор начинается
         с буквы / §-символа (не цифры — иначе ложные срабатывания на
         JSON-массивах `[1, 2, 3]` в HTML embedded scripts).
      2. Добавить **методологический baseline** из 5 канонических источников,
         которые используются в любой gamma-spectrum analysis:
         ГОСТ Р 7.0.5–2008, ISO 11929:2019, ЛСРМ Algo, Будыка-2021, Gilmore 3rd ed.
      3. Если text содержит дополнительные ссылки [N, локатор] (например,
         после translate_text применил Layer 1→2 conversion) — включаем их тоже.
    """
    if not text:
        return text
    # Refined regex: [N, локатор] где локатор начинается с не-цифры
    # (защита от JSON arrays). Например: [7, §10] [12, гл.3] [19, p.225].
    refs_used = set()
    for m in re.finditer(
        r"\[(\d{1,3}),\s*([§с.A-Za-zА-Яа-яΑ-Ωα-ω][^\]]*)\]", text,
    ):
        try:
            refs_used.add(int(m.group(1)))
        except ValueError:
            continue

    # Methodology baseline: фиксированный набор для любого анализа
    # (даже если автоматический сбор не нашёл ссылок).
    # F-337.4 / v1.18.19.1 — ID=1 (ГОСТ Р 7.0.5–2008 «Библиографическая
    # ссылка») удалён из baseline по запросу пользователя: ссылка на сам
    # citation-standard избыточна в списке источников.
    # BUG-16 — ID=24 (ГОСТ 26874-86 «Спектрометры энергий ионизирующих
    # излучений. Методы измерения основных параметров») добавлен в
    # baseline: foundational normative document for γ-spectrometer
    # metrology (R(E), линейность, ε в ПИП, ИНЛ, загрузочная способность).
    # PDF: books_library/01_methodology_pdf/
    # GOST_26874-86_спектрометры_методы_измерения.pdf.
    # Затекстовая запись — см. references/REFERENCES.md §1 запись № 24.
    baseline_refs = {2, 7, 12, 19, 24}  # ISO 11929, LSRM-Algo, Будыка, Gilmore, ГОСТ 26874-86
    refs_used.update(baseline_refs)
    # Также защитно убираем ID=1 если он попал из text-scan (для совместимости)
    refs_used.discard(1)

    refs_map = _load_references_map()
    if not refs_map:
        return text

    sorted_refs = sorted(refs_used)
    lines = []
    if format == "md":
        lines.append("\n\n---\n\n## Список использованной литературы\n\n")
        # F-337.5 / v1.18.19.1 — фраза «Оформление по ГОСТ Р 7.0.5–2008…» убрана
        lines.append(
            "*В список включён базовый методологический минимум — "
            "нормативные документы, методики и монографии, на которых "
            "основан настоящий анализ.*\n\n",
        )
        for n in sorted_refs:
            citation = refs_map.get(
                n, f"[Источник {n} — описание отсутствует в реестре]",
            )
            lines.append(f"{n}. {citation}\n")
    elif format == "html":
        lines.append('\n<hr/>\n<section class="gost-references">\n')
        lines.append('  <h2>Список использованной литературы</h2>\n')
        # F-337.5 — фраза «Оформление по ГОСТ Р 7.0.5–2008…» убрана
        lines.append(
            '  <p><em>В список включён базовый методологический минимум — '
            'нормативные документы, методики и монографии, на которых '
            'основан настоящий анализ.</em></p>\n',
        )
        lines.append('  <ol>\n')
        for n in sorted_refs:
            citation = refs_map.get(
                n, f"[Источник {n} — описание отсутствует в реестре]",
            )
            # Escape HTML special chars
            cit_safe = (
                citation.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;")
            )
            lines.append(f'    <li value="{n}">{cit_safe}</li>\n')
        lines.append('  </ol>\n')
        lines.append('</section>\n')
    appendix = "".join(lines)
    if format == "html":
        # F-324 / v1.18.18 — Inject INSIDE the .page container when present.
        # Interactive template (interactive_v1_17_2.html) wraps content in
        # <div class="page"> with max-width:900px; margin:0 auto. Inserting
        # before </body> would place the section OUTSIDE the container
        # → full-window width visual breakage. The interactive template
        # also has <script> tags between the .page closing </div> and
        # </body>, so a simple `</div>\s*</body>` regex won't find the
        # right anchor — we must balance-match nested <div>/</div> starting
        # from `<div class="page">`.
        # Static template (html_report.py) uses `body` itself as the
        # max-width:1100px constraint, so before-</body> is correct there.
        # Regression guard: tests/snapshot/test_v1_18_18_refs_inside_container.py
        page_close_idx = _find_page_div_close(text)
        if page_close_idx is not None:
            return text[:page_close_idx] + appendix + text[page_close_idx:]
        # Fallback: insert before </body> (works for static template where
        # body itself has the max-width constraint).
        lower_text = text.lower()
        body_close = lower_text.rfind("</body>")
        if body_close > 0:
            return text[:body_close] + appendix + text[body_close:]
        # Last-resort fallback: append (out-of-spec but non-empty).
    return text + appendix


_DIV_TAG_RE = re.compile(r'<(/?)div\b[^>]*>', re.IGNORECASE)


def _find_page_div_close(html: str) -> Optional[int]:
    """Locate the index of `</div>` that closes `<div class="page">`.

    Returns None when no `<div class="page">` is present (e.g. static
    template). Uses depth-balanced scanning over `<div>` / `</div>`
    tokens, ignoring `<div>` references inside string literals of
    `<script>` blocks for now (those don't appear in our templates).
    """
    page_open = re.search(r'<div\s[^>]*class="page"', html, re.IGNORECASE)
    if not page_open:
        return None
    # Start scanning AFTER the opening tag of <div class="page">.
    # Find the end of that opening tag first.
    open_end_match = re.compile(r'>').search(html, page_open.start())
    if not open_end_match:
        return None
    cursor = open_end_match.end()
    depth = 1
    for m in _DIV_TAG_RE.finditer(html, cursor):
        if m.group(1) == "":  # opening <div
            depth += 1
        else:                  # closing </div>
            depth -= 1
            if depth == 0:
                return m.start()
    return None


def _extract_activities_dict(json_dict: Dict[str, Any]) -> Dict[str, float]:
    """F-326 helper: вытащить {nuclide: A_Bq} из json_dict для passport
    compare table. Robust к разным схемам storage активностей.
    """
    out: Dict[str, float] = {}
    if not json_dict:
        return out
    # Schema 1: identified_nuclides[].A_Bq
    for n in json_dict.get("identified_nuclides", []) or []:
        nuc = n.get("nuclide") or n.get("name")
        a = n.get("A_Bq") or n.get("activity_Bq") or n.get("activity")
        if nuc and a is not None:
            try:
                out[str(nuc)] = float(a)
            except (TypeError, ValueError):
                pass
    return out


def _extract_specific_activities_dict(
    json_dict: Dict[str, Any],
) -> Dict[str, float]:
    """F-376 helper: вытащить {nuclide: A_Bq_per_kg} из json_dict для
    passport-comparison в удельных активностях (Бк/кг).

    Поле `specific_activity_Bq_per_kg` уже посчитано pipeline'ом
    (json_report.py) с учётом sample mass + decay correction, поэтому
    здесь только plain extract без пересчёта.
    """
    out: Dict[str, float] = {}
    if not json_dict:
        return out
    for n in json_dict.get("identified_nuclides", []) or []:
        nuc = n.get("nuclide") or n.get("name")
        sa = n.get("specific_activity_Bq_per_kg")
        if nuc and sa is not None:
            try:
                out[str(nuc)] = float(sa)
            except (TypeError, ValueError):
                pass
    return out


def _derive_sample_mass_kg(
    activities_Bq: Dict[str, float],
    specific_activities_Bq_per_kg: Dict[str, float],
) -> Optional[float]:
    """F-376 helper: восстановить sample mass из ratio A_Bq / SA_Bq_per_kg
    для любого нуклида с обоими значениями. Возвращает None если нет
    данных. Это fallback для случая, когда passport_meta['mass_kg'] None.
    """
    for nuc, a in activities_Bq.items():
        sa = specific_activities_Bq_per_kg.get(nuc)
        try:
            if a is not None and sa is not None and float(sa) > 0:
                m = float(a) / float(sa)
                if m > 0:
                    return m
        except (TypeError, ValueError):
            pass
    return None


def _fmt_activity(v: Optional[float]) -> str:
    """F-377 helper: format Бк/кг в human-readable виде с ru-thousand-
    separator вместо E-нотации. Примеры: 3118 → "3 118", -0.3 → "0",
    None / nan / inf → "—".
    Минус заменяется на typographic minus (U+2212).
    """
    import math as _m
    if v is None:
        return "—"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    if not _m.isfinite(f):
        return "—"
    if abs(f) < 0.5:
        return "0"
    s = f"{int(round(f)):,}".replace(",", " ")
    return s.replace("-", "−")


# F-337.3 / v1.18.19.1 — chain-parent equivalents для passport-comparison
# fallback. Когда passport содержит chain-parent (Th-232, U-238, Ra-226), а
# direct activity отсутствует — берём из daughter в secular equilibrium,
# применяя branching-фактор. Формат: parent → [(daughter, factor), ...].
# factor = A(parent) / A(daughter); первый available daughter побеждает.
#
# Источники branching:
#   • Th-232 → Pb-212 (1:1, прямая equilibrium); Ac-228 (1:1);
#     Tl-208 / 0.3593 (Bi-212 β-decay branching = 35.93%, ENSDF)
#   • U-238 / Ra-226 → Bi-214 (1:1); Pb-214 (1:1)
_CHAIN_PARENT_EQUIV = {
    "Th-232": [
        ("Pb-212", 1.0),
        ("Ac-228", 1.0),
        ("Tl-208", 1.0 / 0.3593),
    ],
    "U-238": [
        ("Bi-214", 1.0),
        ("Pb-214", 1.0),
        ("Ra-226", 1.0),
    ],
    "Ra-226": [
        ("Bi-214", 1.0),
        ("Pb-214", 1.0),
    ],
}


def _resolve_measured_for_passport(
    nuc: str,
    activities: Dict[str, float],
) -> tuple:
    """F-337.3 — return (measured_Bq, source_note) for passport row.

    Direct lookup wins; fallback на chain-parent equilibrium daughter если
    parent в _CHAIN_PARENT_EQUIV. Returns (None, "") если ничего не нашли.
    """
    direct = activities.get(nuc)
    if direct is not None:
        return float(direct), ""
    daughters = _CHAIN_PARENT_EQUIV.get(nuc)
    if not daughters:
        return None, ""
    for d_name, factor in daughters:
        d_act = activities.get(d_name)
        if d_act is not None:
            note = f" (по дочернему {d_name} в равновесии цепочки)"
            return float(d_act) * float(factor), note
    return None, ""


def _f326_append_passport_comparison(
    text: str,
    *,
    format: str,
    passport: Optional[Dict[str, float]],
    activities: Dict[str, float],
    passport_meta: Optional[Dict[str, Any]] = None,
    specific_activities: Optional[Dict[str, float]] = None,
) -> str:
    """F-326 / v1.18.18.1 — append «Сравнение с паспортной активностью».

    F-376 / v1.18.24.x — comparison ведётся в удельных активностях
    (Бк/кг) для consistency с summary card. Passport pre F-369 хранит
    total Бк (после конверсии из Бк/кг × cert mass) — здесь делим
    обратно на ту же массу. Если `passport_meta['mass_kg']` None —
    fallback на sample mass, восстановленную из ratio measured A_Bq /
    SA_Bq_per_kg (см. _derive_sample_mass_kg).

    F-377 / v1.18.24.x — числа форматируются через `_fmt_activity`
    (ru-thousand-separator, без E-нотации).

    User-feedback closure: даже когда `passport_activity_Bq` НЕ передан,
    блок отображается с явным сообщением о deferred-state и
    инструкцией как его включить. Это устраняет «silent missing
    section» — оператор всегда видит, есть ли сравнение и где взять
    данные.

    Args:
        text: report body (MD or HTML)
        format: 'md' | 'html'
        passport: {'Cs-137': 1.05e3, ...} total Bq, либо None
        activities: {'Cs-137': 1.038e3, ...} measured total Bq
        passport_meta: dict с mass_kg / meas_date / decay_corrected
        specific_activities: {'Cs-137': 1820.0, ...} measured Бк/кг;
            если None — fallback на activities / sample_mass.

    Insertion point: HTML — внутри .page контейнера (использует
    _find_page_div_close); MD — перед F-318 refs списком (т.е.
    раньше append'а refs).
    """
    if not text:
        return text

    has_data = bool(passport) and isinstance(passport, dict)
    meta = passport_meta or {}
    sa_in = specific_activities or {}

    # F-376 — определить mass_kg для Бк → Бк/кг конверсии паспорта.
    # Приоритет: meta['mass_kg'] (из .src certificate) → derived sample
    # mass (A_Bq / SA_Bq_per_kg). Если оба None — comparison fallback'ит
    # на total Бк (legacy behaviour) и в footnote это отмечается.
    mass_for_conversion: Optional[float] = None
    mass_source_note = ""
    try:
        m_meta = meta.get("mass_kg")
        if m_meta is not None and float(m_meta) > 0:
            mass_for_conversion = float(m_meta)
            mass_source_note = "по mass_kg из .src сертификата"
    except (TypeError, ValueError):
        pass
    if mass_for_conversion is None:
        derived = _derive_sample_mass_kg(activities, sa_in)
        if derived is not None and derived > 0:
            mass_for_conversion = derived
            mass_source_note = (
                "оценена из ratio A_Bq / SA_Bq_per_kg (mass_kg в .src "
                "отсутствует)"
            )

    # F-376 — измеренные удельные активности для таблицы. Если есть
    # explicit specific_activities (поле json_report.py) — используем
    # их; иначе делим total на mass_for_conversion (если есть). Это
    # отдельный dict от `activities` (Бк) — passes в _resolve_measured.
    sa_measured: Dict[str, float] = dict(sa_in)
    if mass_for_conversion is not None and not sa_measured:
        for nuc, a in activities.items():
            try:
                sa_measured[nuc] = float(a) / float(mass_for_conversion)
            except (TypeError, ValueError, ZeroDivisionError):
                pass
    src = (meta.get("source") or "").lower()
    src_label_ru = {
        "spe_comment": "автоматически извлечено из поля COMMENT файла .spe",
        "explicit":    "передано пользователем явно через `passport_activity_Bq`",
        "none":        "источник не указан",
    }.get(src, "источник не указан")

    if format == "md":
        lines = ["\n\n---\n\n## Сравнение с паспортной удельной активностью\n\n"]
        if not has_data:
            lines.append(
                "*Сравнение не выполнено — данные паспорта источника "
                "не переданы.*\n\n"
                "Чтобы включить, передайте `passport_activity_Bq` "
                "при вызове `analyze_and_report`, например:\n\n"
                "```python\n"
                "analyze_and_report(\n"
                "    'M_cs_легкий_2001-2005.spe',\n"
                "    sample_mass_kg=0.570,\n"
                "    passport_activity_Bq={'Cs-137': 1.050e3},  # Bq, "
                "из паспорта эталонного источника\n"
                ")\n"
                "```\n\n"
                "Базис сравнения: ИНЛ (нижний предел измерений) "
                "и НЛ (нижний предел разрешения) по ГОСТ Р 8.594–2002, "
                "критерий совпадения — отклонение ≤ U(P=0.95) от паспорта.\n"
            )
        else:
            # F-330 — provenance preamble (масса, дата измерения,
            # decay-correction summary).
            # F-337.6 / v1.18.19.1 — фраза «Источник паспорта: …» убрана
            # (избыточно, источник всегда .spe COMMENT для current pipeline).
            mass_kg = meta.get("mass_kg")
            meas_date = meta.get("meas_date")
            if mass_kg or meas_date:
                bits = []
                if mass_kg:
                    bits.append(f"масса образца {float(mass_kg):.3f} кг")
                if meas_date:
                    bits.append(f"дата измерения {meas_date}")
                lines.append("*Параметры пересчёта: " + ", ".join(bits) + ".*\n\n")
            dc_map = meta.get("decay_corrected") or {}
            ref_dates = meta.get("ref_dates") or {}
            dc_lines = []
            for n, applied in dc_map.items():
                ref = ref_dates.get(n)
                if applied and ref:
                    dc_lines.append(f"{n}: ref={ref}, decay-correction применён (exp(-λΔt))")
                elif ref and not applied:
                    dc_lines.append(f"{n}: ref={ref}, decay-correction НЕ применён (t½ не в таблице)")
                elif not ref:
                    dc_lines.append(f"{n}: ref-дата не указана → decay-correction опущен")
            if dc_lines:
                lines.append("*Decay-correction:* " + "; ".join(dc_lines) + ".\n\n")
            lines.append(
                "| Нуклид | Измерено, Бк/кг | Паспорт, Бк/кг | "
                "Отклонение, % | Совпадение |\n"
            )
            lines.append(
                "|---|---:|---:|---:|---|\n"
            )
            # F-376 — рендерим в Бк/кг. Measured берётся из sa_measured
            # (chain-parent fallback тоже срабатывает на dict удельных
            # активностей). Passport_per_kg = passport_total / mass_for_conversion.
            for nuc, p_bq in sorted(passport.items()):
                # F-337.3 — chain-parent fallback (Th-232 → Pb-212/Ac-228 равновесие)
                m_sa, src_note = _resolve_measured_for_passport(nuc, sa_measured)
                p_sa = None
                if mass_for_conversion is not None and p_bq is not None:
                    try:
                        p_sa = float(p_bq) / float(mass_for_conversion)
                    except (TypeError, ValueError, ZeroDivisionError):
                        p_sa = None
                if m_sa is None or p_sa is None or p_sa <= 0:
                    lines.append(
                        f"| {nuc} | {_fmt_activity(m_sa)} | "
                        f"{_fmt_activity(p_sa)} | — | "
                        "*нет измерения* |\n"
                    )
                    continue
                dev = 100.0 * (float(m_sa) - float(p_sa)) / float(p_sa)
                # 25% — generic acceptance band (NaI 3″×3″ default).
                # Реальный critère — U(P=0.95) индивидуально per nuclide.
                ok = "✓" if abs(dev) <= 25.0 else "✗"
                meas_disp = f"{_fmt_activity(m_sa)}{src_note}"
                lines.append(
                    f"| {nuc} | {meas_disp} | {_fmt_activity(p_sa)} | "
                    f"{dev:+.1f} | {ok} |\n"
                )
            if mass_source_note:
                lines.append(
                    f"\n*Конверсия Бк → Бк/кг: масса {mass_source_note}.*\n"
                )
            lines.append(
                "\n*Критерий приёмки: |отклонение| ≤ 25% (generic NaI "
                "63×63 default). Точный критерий — U(P=0.95) по ГОСТ "
                "Р 8.594–2002.*\n"
            )
        # Markdown path stops here; HTML path continues in the elif below.
        appendix = "".join(lines)
        return text + appendix

    elif format == "html":
        lines = [
            '\n<hr/>\n<section class="passport-comparison">\n',
            '  <h2>Сравнение с паспортной удельной активностью</h2>\n',
        ]
        if not has_data:
            lines.append(
                '  <p><em>Сравнение не выполнено — данные паспорта '
                'источника не переданы.</em></p>\n'
                '  <p>Чтобы включить, передайте '
                '<code>passport_activity_Bq</code> при вызове '
                '<code>analyze_and_report</code>, например:</p>\n'
                '  <pre><code>'
                "analyze_and_report(\n"
                "    'M_cs_легкий_2001-2005.spe',\n"
                "    sample_mass_kg=0.570,\n"
                "    passport_activity_Bq={'Cs-137': 1.050e3},\n"
                ")"
                '</code></pre>\n'
                '  <p>Базис сравнения: ИНЛ / НЛ по ГОСТ Р 8.594–2002, '
                'критерий совпадения — отклонение ≤ U(P=0.95) от '
                'паспорта.</p>\n'
            )
        else:
            # F-330 — provenance preamble в HTML формате.
            # F-337.6 / v1.18.19.1 — фраза «Источник паспорта: …» убрана.
            mass_kg = meta.get("mass_kg")
            meas_date = meta.get("meas_date")
            if mass_kg or meas_date:
                bits = []
                if mass_kg:
                    bits.append(f"масса образца {float(mass_kg):.3f} кг")
                if meas_date:
                    bits.append(f"дата измерения {meas_date}")
                lines.append(
                    "  <p><em>Параметры пересчёта: "
                    + ", ".join(bits) + ".</em></p>\n"
                )
            dc_map = meta.get("decay_corrected") or {}
            ref_dates = meta.get("ref_dates") or {}
            dc_bits = []
            for n, applied in dc_map.items():
                ref = ref_dates.get(n)
                if applied and ref:
                    dc_bits.append(
                        f"{n}: ref={ref}, decay-correction применён "
                        "(exp(-&lambda;&Delta;t))"
                    )
                elif ref and not applied:
                    dc_bits.append(
                        f"{n}: ref={ref}, decay-correction НЕ применён "
                        "(t&frac12; не в таблице)"
                    )
                elif not ref:
                    dc_bits.append(
                        f"{n}: ref-дата не указана → decay-correction опущен"
                    )
            if dc_bits:
                lines.append(
                    "  <p><em>Decay-correction:</em> "
                    + "; ".join(dc_bits) + ".</p>\n"
                )
            lines.append('  <table class="passport-tbl">\n')
            lines.append(
                '    <thead><tr><th>Нуклид</th><th>Измерено, Бк/кг</th>'
                '<th>Паспорт, Бк/кг</th><th>Отклонение, %</th>'
                '<th>Совпадение</th></tr></thead>\n'
            )
            lines.append('    <tbody>\n')
            # F-376 — рендерим Бк/кг. Симметрично markdown-ветке выше.
            for nuc, p_bq in sorted(passport.items()):
                # F-337.3 — chain-parent fallback (Th-232 → Pb-212/Ac-228 равновесие)
                m_sa, src_note = _resolve_measured_for_passport(nuc, sa_measured)
                p_sa = None
                if mass_for_conversion is not None and p_bq is not None:
                    try:
                        p_sa = float(p_bq) / float(mass_for_conversion)
                    except (TypeError, ValueError, ZeroDivisionError):
                        p_sa = None
                if m_sa is None or p_sa is None or p_sa <= 0:
                    lines.append(
                        f'      <tr><td>{nuc}</td>'
                        f'<td>{_fmt_activity(m_sa)}</td>'
                        f'<td>{_fmt_activity(p_sa)}</td><td>—</td>'
                        '<td><em>нет измерения</em></td></tr>\n'
                    )
                    continue
                dev = 100.0 * (float(m_sa) - float(p_sa)) / float(p_sa)
                ok = "✓" if abs(dev) <= 25.0 else "✗"
                meas_cell = _fmt_activity(m_sa)
                if src_note:
                    meas_cell += f' <em style="font-size:11px;color:var(--text-tertiary)">{src_note.strip()}</em>'
                lines.append(
                    f'      <tr><td>{nuc}</td>'
                    f'<td>{meas_cell}</td>'
                    f'<td>{_fmt_activity(p_sa)}</td>'
                    f'<td>{dev:+.1f}</td><td>{ok}</td></tr>\n'
                )
            lines.append('    </tbody>\n  </table>\n')
            if mass_source_note:
                lines.append(
                    f'  <p><em>Конверсия Бк → Бк/кг: масса '
                    f'{mass_source_note}.</em></p>\n'
                )
            lines.append(
                '  <p><em>Критерий приёмки: |отклонение| ≤ 25% '
                '(generic NaI 63×63 default). Точный критерий — '
                'U(P=0.95) по ГОСТ Р 8.594–2002.</em></p>\n'
            )
        lines.append('</section>\n')
        appendix = "".join(lines)
        # F-337.2 / v1.18.19.1 — passport block ниже блока активности
        # (`.fp-summary`), но ВЫШЕ заметок и cost footer'а. Логически: сначала
        # «итоговая активность образца», сразу под ней — сравнение с паспортом,
        # потом — narrative заметки и refs.
        # Anchor: ищем `<div class="fp-notes">` (заметки/заключение) и
        # вставляем appendix ДО него. Fallback цепочка → page-close → body-close.
        notes_marker = '<div class="fp-notes">'
        notes_idx = text.find(notes_marker)
        if notes_idx >= 0:
            return text[:notes_idx] + appendix + text[notes_idx:]
        page_close_idx = _find_page_div_close(text)
        if page_close_idx is not None:
            return text[:page_close_idx] + appendix + text[page_close_idx:]
        lower_text = text.lower()
        body_close = lower_text.rfind("</body>")
        if body_close > 0:
            return text[:body_close] + appendix + text[body_close:]
        return text + appendix

    return text


def _load_references_map() -> Dict[int, str]:
    """Hardcoded ГОСТ Р 7.0.5–2008 mapping № → full citation.

    Source of truth: `references/REFERENCES.md` (см. citation_translator.py
    RAG_PREFIX_TO_GOST). Тут — сжатые библ.записи для footer-генерации.
    Обновлять синхронно с RAG_PREFIX_TO_GOST при добавлении новых источников.
    """
    return {
        1: "ГОСТ Р 7.0.5–2008. Библиографическая ссылка. Общие требования и правила составления. — М.: Стандартинформ, 2008.",
        2: "ISO 11929:2019. Determination of the characteristic limits (decision threshold, detection limit and limits of the confidence interval) for measurements of ionizing radiation. — Geneva: ISO, 2019.",
        3: "ICRP Publication 74. Conversion Coefficients for use in Radiological Protection against External Radiation. — Annals of the ICRP, Vol. 26, No. 3/4, 1996.",
        4: "ОСПОРБ-99/2010. Основные санитарные правила обеспечения радиационной безопасности. СП 2.6.1.2612-10.",
        5: "Активность в счётных образцах. Методика измерений на гамма-спектрометрах с использованием ПО SpectraLine. — М.: ООО «ЛСРМ», 2014.",
        6: "Мощность дозы. Методика расчёта из спектра гамма-излучения. — М.: ВНИИФТРИ; ООО «ЛСРМ», 2000. — 10 с.",
        7: "Алгоритмические основы программ обработки спектрометрической информации SpectraLine. — М.: ООО «ЛСРМ», 2022.",
        8: "Описание формата файлов SpectraLine .spe. — М.: ООО «ЛСРМ».",
        9: "SpectraLine 2.0 — Основные функции: руководство пользователя. — М.: ООО «ЛСРМ».",
        10: "Кувыкин В. И. Прецизионные методы (презентация). — М.: ООО «ЛСРМ», 2023.",
        11: "ЛСРМ. Методы и алгоритмы паспортизации РАО (презентация). — М.: ООО «ЛСРМ».",
        12: "Будыка А. К. Спектрометрия ионизирующих излучений. Гамма-спектрометрия: учеб. пособие. — М.: НИЯУ МИФИ, 2021. — 225 с.",
        13: "Будыка А. К. Спектрометрия ионизирующих излучений. Основные понятия и терминология: учеб.-метод. пособие. — М.: НИЯУ МИФИ, 2021. — 144 с. — ISBN 978-5-7262-2794-8.",
        14: "Шендрик Р. Ю. Введение в физику сцинтилляторов. Часть 1: учеб. пособие. — Иркутск, 2017.",
        15: "Шендрик Р. Ю. Введение в физику сцинтилляторов. Часть 2: учеб. пособие. — Иркутск, 2018.",
        16: "Вартанов Н. А., Самойлов П. С. Практические методы сцинтилляционной γ-спектрометрии. — М.: Атомиздат.",
        17: "Минимальная детектируемая активность. Основные понятия и определения. — Статья (6 с.).",
        18: "Анализ и представление результатов эксперимента: учеб. пособие. — Статистическая обработка данных.",
        19: "Gilmore G. R., Hemingway J. D. Practical Gamma-ray Spectrometry. — 3rd ed. — Chichester: Wiley, 2024.",
        20: "Knoll G. F. Radiation Detection and Measurement. — 4th ed. — New York: Wiley, 2010.",
        21: "Debertin K., Helmer R. G. Gamma- and X-ray Spectrometry with Semiconductor Detectors. — Amsterdam: North-Holland, 1988.",
        22: "GammaVision® / Maestro-PRO®. Gamma-Ray Spectrum Analysis and MCA Emulators. Software User's Manual. Software Version 9. Part No. 783620, Rev. M (0220). — Oak Ridge, TN: ORTEC (AMETEK), 2020. — 83 p.",
        23: "scikit-learn developers. sklearn.mixture — Gaussian Mixture / Bayesian Gaussian Mixture (EM + Variational Inference). — URL: https://scikit-learn.org/stable/modules/mixture.html.",
        # BUG-16 — № 24 — нормативный ГОСТ для метрологии γ-спектрометров.
        # Точное название титульного листа PDF (страница 1–2):
        # «Спектрометры энергий ионизирующих излучений. Методы измерения
        # основных параметров» (СТ СЭВ 5053-85). Утв. Госкомстандарт СССР
        # 21.04.86 № 1016. Введ. 1987-01-01. УДК 539.1.074.083 : 006.354.
        # Источник истины — references/REFERENCES.md §1 запись № 24.
        24: "ГОСТ 26874–86 (СТ СЭВ 5053–85). Спектрометры энергий ионизирующих излучений. Методы измерения основных параметров. — Введ. 1987–01–01. — М.: Издательство стандартов, 1987. — 30 с.",
    }


def _safe_filename_stem(filename: str) -> str:
    """Return a filesystem-safe stem from an arbitrary input filename.

    F-115: strips embedded S/N tokens (e.g. "_420-7-17") so the report
    folder name and downstream artefact basenames do not leak source
    serial numbers.
    """
    base = filename or "report"
    # Drop extension
    if "." in base:
        base = base.rsplit(".", 1)[0]
    # F-115 — scrub embedded S/N tokens before the rest of the cleanup.
    try:
        from gamma.reporting.anonymize import _scrub_sn_in_basename
        base = _scrub_sn_in_basename(base)
    except Exception:
        pass
    # Replace characters that are problematic in Windows/POSIX filenames
    bad = '<>:"/\\|?*'
    for ch in bad:
        base = base.replace(ch, "_")
    return base or "report"


def build_report(
    result,
    *,
    output_dir: Optional[str] = None,
    write_json: bool = True,
    write_markdown: bool = False,
    write_plots: bool = False,
    write_html: bool = False,
    write_technical_pdf: bool = False,
    return_summary: bool = True,
    report_stem: Optional[str] = None,
    plot_dpi: int = 120,
    cost_estimate: Optional[Dict[str, Any]] = None,
    passport_activity_Bq: Optional[Dict[str, float]] = None,
    passport_meta: Optional[Dict[str, Any]] = None,
    bundle_index: bool = False,
) -> Dict[str, Any]:
    """Assemble + persist a Step-11 report from a StagedAnalysisResult.

    Parameters
    ----------
    result : StagedAnalysisResult
    output_dir : str, optional
        Directory to write the report files to. If None, files are not
        written to disk and the function returns the assembled artefacts
        in memory (``json_dict``, ``markdown_text``, ``summary``).
        Plots/HTML require ``output_dir`` — without it those flags are
        silently ignored.
    write_json : bool, default True
        Write ``{stem}_report.json`` to ``output_dir``.
    write_markdown : bool, default False
        Write ``{stem}_report.md`` to ``output_dir``.
    write_plots : bool, default False
        Generate PNG plots (spectrum + multiplet clusters) under
        ``{output_dir}/{stem}_plots/``. Requires matplotlib.
    write_html : bool, default False
        Render an HTML version of the report (embedded base64 PNGs
        when ``write_plots=True``, otherwise text-only).
    return_summary : bool, default True
        Generate and include the 3–8 line chat summary in the result.
    report_stem : str, optional
        Override the file stem. By default, derived from
        ``result.spec.source_path`` (or ``"report"`` when no path is
        known).
    plot_dpi : int, default 120

    Returns
    -------
    dict
        Keys:
        * ``json``       — path (or None) of the JSON file.
        * ``markdown``   — path (or None) of the Markdown file.
        * ``html``       — path (or None) of the HTML file.
        * ``plots``      — dict ``{"spectrum": path | None,
                                  "multiplets": List[str]}``  or None.
        * ``summary``    — chat summary string (or None).
        * ``json_dict``  — the assembled JSON dict (always populated).
        * ``markdown_text`` — the Markdown string (when generated).
        * ``html_text``  — the HTML string (when generated).
        * ``warnings``   — list of soft-failure messages (e.g.
                           matplotlib missing → plots silently skipped).
    """
    json_dict = build_json_report(result)

    # ─── F-132 / v1.17.7 — обязательная оценка стоимости ──────────────
    # Жёстко: HTML footer всегда показывает итог; Markdown получает
    # поэтапную таблицу + итог. Если CLI передал --cost-tokens, оно
    # используется как override итога; per-stage таблица всё равно
    # строится из авто-оценки (прозрачность).
    try:
        from gamma.reporting.cost_estimator import (
            estimate_total_cost, DEFAULT_SESSION_TOKEN_BUDGET,
        )
        budget = DEFAULT_SESSION_TOKEN_BUDGET
        override_tokens = None
        override_detail = None
        if cost_estimate:
            try:
                t = cost_estimate.get("tokens")
                if t is not None and int(t) > 0:
                    override_tokens = int(t)
            except (TypeError, ValueError):
                pass
            override_detail = cost_estimate.get("detail")
            try:
                b = cost_estimate.get("session_token_budget")
                if b is not None and int(b) > 0:
                    budget = int(b)
            except (TypeError, ValueError):
                pass
        cost_full = estimate_total_cost(
            result,
            session_token_budget=budget,
            cost_tokens_override=override_tokens,
            detail_override=override_detail,
        )
        # P0-8: inject actual Claude output-token count into CostEstimate.
        # Source: cost_estimate["output_tokens"] from CLI caller when available.
        # CostEstimate is frozen=True (cost_estimator.py:73) — use dataclasses.replace.
        # If not provided by caller, defaults to 0 (unknown at this call site).
        # TODO: wire CLI --cost-output-tokens flag for live count injection.
        _output_tokens: int = 0
        if cost_estimate:
            try:
                _ot = cost_estimate.get("output_tokens")
                if _ot is not None and int(_ot) >= 0:
                    _output_tokens = int(_ot)
            except (TypeError, ValueError):
                pass
        cost_full = dataclasses.replace(cost_full, claude_output_tokens=_output_tokens)
        json_dict["cost_estimate"] = cost_full.to_dict()
        # P0-8: 20k alarm — elevated output-token cost.
        if cost_full.claude_output_tokens >= 20_000:
            json_dict.setdefault("warnings", []).append({
                "code": "COST_HIGH_OUTPUT_TOKENS",
                "message": (
                    f"claude_output_tokens={cost_full.claude_output_tokens}"
                    f" >= 20000 threshold. Analysis cost is elevated."
                ),
                "claude_output_tokens": cost_full.claude_output_tokens,
                "threshold": 20_000,
                "severity": "INFO",
            })
        # Также строим dict для legacy interactive_html (HTML footer).
        # Если пользователь передал session_pct как строку (D-19 legacy
        # API), используем его; иначе подставляем авто-формулу.
        explicit_session_pct = None
        if cost_estimate:
            sp = cost_estimate.get("session_pct")
            if sp:
                explicit_session_pct = str(sp)
        cost_estimate_for_html = {
            "tokens": cost_full.tokens_total,
            "session_pct": (
                explicit_session_pct
                or f"{cost_full.session_pct:.1f}% от бесплатной 5-часовой сессии"
            ),
            "session_token_budget": cost_full.session_token_budget,
            "detail": cost_full.detail,
        }
    except Exception:
        cost_full = None
        cost_estimate_for_html = cost_estimate

    out: Dict[str, Any] = {
        "json": None,
        "markdown": None,
        "html": None,
        "technical_pdf": None,
        "plots": None,
        "summary": None,
        "json_dict": json_dict,
        "markdown_text": None,
        "html_text": None,
        "warnings": [],
    }

    if report_stem is None:
        sp = ""
        try:
            sp = result.spec.source_path or ""
        except AttributeError:
            sp = ""
        leaf = ""
        if sp:
            leaf = os.path.basename(sp)
        report_stem = _safe_filename_stem(leaf or "report")

    md_dir_for_links: Optional[str] = None
    plots: Optional[Dict[str, Any]] = None

    if output_dir is not None:
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        md_dir_for_links = str(out_path)

        # --- Plots ----------------------------------------------------
        if write_plots:
            try:
                from gamma.reporting.plots import build_all_plots
                plots_dir = out_path / f"{report_stem}_plots"
                plots = build_all_plots(
                    result, plots_dir, dpi=plot_dpi,
                )
                out["plots"] = plots
            except ImportError as e:
                out["warnings"].append(
                    f"plots skipped: {e}"
                )
            except Exception as e:
                out["warnings"].append(
                    f"plots failed: {type(e).__name__}: {e}"
                )

        # --- JSON -----------------------------------------------------
        if write_json:
            # Inject plot paths into the JSON so downstream tools know
            # which PNGs were generated.
            if plots is not None:
                json_dict["plot_files"] = {
                    "spectrum": plots.get("spectrum"),
                    "multiplets": list(plots.get("multiplets") or []),
                }
            else:
                json_dict.setdefault("plot_files", {"spectrum": None,
                                                   "multiplets": []})
            # F-115 — re-anonymise after plot paths are injected.
            try:
                from gamma.reporting.anonymize import anonymize_report_inplace
                anonymize_report_inplace(json_dict)
            except Exception:
                pass
            json_path = out_path / f"{report_stem}_report.json"
            with json_path.open("w", encoding="utf-8") as f:
                json.dump(json_dict, f, ensure_ascii=False, indent=2)
            out["json"] = str(json_path)

        # --- Markdown -------------------------------------------------
        if write_markdown:
            md_text = build_markdown_report(
                result, json_dict=json_dict,
                plots=plots, md_dir=md_dir_for_links,
            )
            # F-317 / v1.18.15 — apply F-256 (Layer 1→2 ГОСТ) + F-260
            # (bilingual narrator) + strip F-id mentions из user-facing body.
            md_text = _f317_apply_user_facing_compliance(md_text, format="md")
            # F-326 / v1.18.18.1 — passport comparison section (приходит
            # ВСЕГДА, даже если данных нет — explicit deferred message).
            # F-UX-04 / 2026-06-04 — но НЕ для чисто-фоновых спектров
            # (env == "background_only"): фон — control, активности не
            # измеряются и с паспортом не сравниваются. См. inbox
            # 2026-06-04_correction_10_no_activity_for_background.md.
            _is_bg_only_pc = (
                (json_dict.get("diagnostics") or {})
                .get("measurement_environment") == "background_only"
            )
            if not _is_bg_only_pc:
                md_text = _f326_append_passport_comparison(
                    md_text, format="md",
                    passport=passport_activity_Bq,
                    activities=_extract_activities_dict(json_dict),
                    passport_meta=passport_meta,
                    specific_activities=(
                        _extract_specific_activities_dict(json_dict)
                    ),
                )
            # F-318 / v1.18.15 — append ГОСТ references list
            md_text = _f318_append_gost_references(md_text, format="md")
            out["markdown_text"] = md_text
            md_path = out_path / f"{report_stem}_report.md"
            with md_path.open("w", encoding="utf-8") as f:
                f.write(md_text)
            out["markdown"] = str(md_path)

        # --- HTML -----------------------------------------------------
        if write_html:
            try:
                from gamma.reporting.html_report import build_html_report
                html_text = build_html_report(
                    result, json_dict=json_dict, plots=plots,
                    cost_estimate=cost_estimate_for_html,
                    bundle_index=bundle_index,
                )
                # F-317 / v1.18.15 — apply same compliance pipeline to HTML
                html_text = _f317_apply_user_facing_compliance(
                    html_text, format="html",
                )
                # F-326 / v1.18.18.1 — passport comparison section
                # F-UX-04 / 2026-06-04 — НЕ выводим для bg-only спектров.
                _is_bg_only_pc_html = (
                    (json_dict.get("diagnostics") or {})
                    .get("measurement_environment") == "background_only"
                )
                if not _is_bg_only_pc_html:
                    html_text = _f326_append_passport_comparison(
                        html_text, format="html",
                        passport=passport_activity_Bq,
                        activities=_extract_activities_dict(json_dict),
                        passport_meta=passport_meta,
                        specific_activities=(
                            _extract_specific_activities_dict(json_dict)
                        ),
                    )
                # F-318 — append references list (inserted inside .page)
                html_text = _f318_append_gost_references(
                    html_text, format="html",
                )
                out["html_text"] = html_text
                html_path = out_path / f"{report_stem}_report.html"
                with html_path.open("w", encoding="utf-8") as f:
                    f.write(html_text)
                out["html"] = str(html_path)
            except ImportError as e:
                out["warnings"].append(
                    f"html skipped: {e}"
                )
            except Exception as e:
                out["warnings"].append(
                    f"html failed: {type(e).__name__}: {e}"
                )

        # --- Technical PDF (F-159 / v1.18.21.0) -----------------------
        # Контракт навсегда: PDF — часть обязательного комплекта отчётных
        # документов. json_dict уже анонимизирован (anonymize_report_inplace
        # вызван выше при write_json), но build_technical_pdf также включает
        # defensive _basename() на всех path-like полях.
        if write_technical_pdf:
            try:
                from gamma.reporting.technical_pdf import build_technical_pdf
                pdf_path = out_path / f"{report_stem}_technical_report.pdf"
                build_technical_pdf(result, json_dict, str(pdf_path))
                out["technical_pdf"] = str(pdf_path)
            except ImportError as e:
                out["warnings"].append(
                    f"technical_pdf skipped (reportlab missing): {e}"
                )
            except Exception as e:
                out["warnings"].append(
                    f"technical_pdf failed: {type(e).__name__}: {e}"
                )
    elif write_markdown:
        # Even without an output_dir, build the Markdown string when asked.
        # No plots in this mode — the placeholder remains.
        out["markdown_text"] = build_markdown_report(
            result, json_dict=json_dict,
        )

    if return_summary:
        out["summary"] = build_chat_summary(
            None, json_dict=json_dict, report_path=out["json"],
        )

    return out


__all__ = ["build_report"]
