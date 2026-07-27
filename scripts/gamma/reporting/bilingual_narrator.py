# -*- coding: utf-8 -*-
"""F-260 — Двуязычный narrative enricher для отчётов (ru/en).

Использует словарь Будыка-2021 (`data/glossary_budyka_2021.json`, 255 терминов,
[13] по ГОСТ Р 7.0.5–2008) для авто-добавления английского эквивалента
в скобках при ПЕРВОМ упоминании русского термина в Markdown/HTML-отчёте.

Контракт (F-260, v1.17.10):
- Любой `reporting/*.py` (markdown_report.py, html_report.py, pdf_export.py),
  генерирующий user-facing narrative, обязан вызвать `enrich_text(content)`
  ПЕРЕД финальной записью.
- JSON-отчёты (json_report.py) НЕ обогащаются (machine-consumed).
- Pipeline: enrich_text(narrative) → translate_text(narrative) (Layer 1→Layer 2)
  → save.

Пример:
    Было: «Обнаружено 12 пиков полного поглощения. Самопоглощение учтено.»
    Стало: «Обнаружено 12 пиков полного поглощения (full-energy peak).
            Самопоглощение (self-absorption) учтено.»

См. также:
- references/GLOSSARY_BUDYKA_2021.md — Markdown глоссарий (для людей)
- data/glossary_budyka_2021.json — машинный JSON (для этого модуля)
- F-108 контракт ru/en narrative
- F-256 двухслойная схема ссылок
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

# Default glossary path относительно репозитория
DEFAULT_GLOSSARY_PATH = Path(__file__).resolve().parents[3] / "data" / "glossary_budyka_2021.json"


_ABBR_RE = re.compile(r"\s*\(([А-ЯA-Z][А-ЯA-Z\d-]{1,15})\)\s*$")


@dataclass(frozen=True)
class GlossaryTerm:
    """Один термин словаря."""
    ru: str
    en: Tuple[str, ...]
    page_from: int
    page_to: int
    id: str = ""

    @property
    def primary_en(self) -> str:
        """Основной английский эквивалент (первый из списка)."""
        return self.en[0] if self.en else ""

    @property
    def match_forms(self) -> Tuple[str, ...]:
        """Все формы которые матчим в тексте.

        Например для «Пик полного поглощения (ППП)»: ('Пик полного поглощения', 'ППП').
        Для «Минимальная детектируемая активность (МДА)»: ('Минимальная детектируемая активность', 'МДА').
        Для одиночного «Активность»: ('Активность',).
        """
        m = _ABBR_RE.search(self.ru)
        if m:
            main = self.ru[:m.start()].strip()
            abbr = m.group(1).strip()
            return (main, abbr) if main else (abbr,)
        return (self.ru,)


def load_glossary(path: Optional[Path] = None) -> List[GlossaryTerm]:
    """Загрузить словарь из JSON.

    Args:
        path: путь к glossary_budyka_2021.json; по умолчанию data/.

    Returns:
        Список GlossaryTerm, отсортированный по длине ru-термина (длинные первыми)
        — чтобы при матчинге сначала находились составные термины («пик полного
        поглощения») и только потом простые («пик»).
    """
    p = path or DEFAULT_GLOSSARY_PATH
    data = json.loads(p.read_text(encoding="utf-8"))
    terms: List[GlossaryTerm] = []
    for t in data.get("terms", []):
        ru = t.get("ru", "").strip()
        if not ru:
            continue
        en = t.get("en", [])
        if isinstance(en, str):
            en = [en]
        terms.append(GlossaryTerm(
            ru=ru,
            en=tuple(s.strip() for s in en if s.strip()),
            page_from=int(t.get("page_from", 0) or 0),
            page_to=int(t.get("page_to", 0) or 0),
            id=t.get("id", ""),
        ))
    # Сортировка по убыванию длины — составные термины раньше простых
    terms.sort(key=lambda x: -len(x.ru))
    return terms


def _is_inside_protected(text: str, start: int) -> bool:
    """True если позиция start внутри code-block / link / уже-обогащённого fragment."""
    before = text[:start]
    # Inside fenced code block? (count ``` before pos)
    if before.count("```") % 2 == 1:
        return True
    # Inside inline code (single backtick)?
    # Считаем количество непарных одиночных backticks
    if before.count("`") - before.count("```") * 3 != 0 and before.rfind("`") > before.rfind("`\n"):
        # heuristic: внутри inline code если последний backtick не закрыт
        last_bt = before.rfind("`")
        nl = before.rfind("\n", 0, last_bt)
        line_after_bt = before[last_bt + 1:]
        if "`" not in line_after_bt:
            return True
    # Inside markdown link [text](url)?
    bracket_depth = 0
    for ch in before[max(0, start - 200):]:
        if ch == "[":
            bracket_depth += 1
        elif ch == "]":
            bracket_depth -= 1
    if bracket_depth > 0:
        return True
    return False


def enrich_text(
    text: str,
    glossary: Optional[List[GlossaryTerm]] = None,
    *,
    max_enrichments_per_term: int = 1,
    min_term_length: int = 4,
    skip_terms: Optional[Iterable[str]] = None,
) -> Tuple[str, Dict[str, int]]:
    """Обогатить текст английскими эквивалентами при первом упоминании.

    Args:
        text: исходный Markdown / plain-text.
        glossary: предзагруженный словарь; если None — load_glossary().
        max_enrichments_per_term: сколько раз обогащать один термин (default 1 — только первое).
        min_term_length: минимальная длина ru-термина для матчинга (фильтрует
            «А», «И» — избежать спама).
        skip_terms: явный skip-list (ru-термины которые НЕ обогащать).

    Returns:
        (enriched_text, statistics) — statistics: {term_ru: count_of_enrichments}.

    Examples:
        >>> g = [GlossaryTerm(ru="самопоглощение", en=("self-absorption",), page_from=96, page_to=96)]
        >>> out, stat = enrich_text("Учтено самопоглощение в матрице.", g)
        >>> out
        'Учтено самопоглощение (self-absorption) в матрице.'
        >>> stat
        {'самопоглощение': 1}
        >>> out, stat = enrich_text("Самопоглощение и снова самопоглощение.", g)
        >>> out
        'Самопоглощение (self-absorption) и снова самопоглощение.'
    """
    if glossary is None:
        glossary = load_glossary()
    skip = set(s.lower() for s in (skip_terms or []))

    counts: Dict[str, int] = {}

    for term in glossary:
        if term.ru.lower() in skip:
            continue
        if not term.primary_en:
            continue
        # Каждая форма (полная + аббревиатура) — собственный счётчик
        for form in term.match_forms:
            if len(form) < min_term_length and not (form.isupper() and len(form) >= 2):
                # Пропускаем короткие, но разрешаем UPPERCASE аббревиатуры (МДА, ППП)
                continue
            pattern = re.compile(
                r"(?<![А-Яа-яA-Za-z\-])" + re.escape(form) + r"(?![А-Яа-яA-Za-z\-])",
                re.IGNORECASE if not form.isupper() else 0,
            )
            n_done = [0]
            text_state = [text]

            def _sub(m: re.Match, n_done=n_done, term=term, text_state=text_state) -> str:
                if n_done[0] >= max_enrichments_per_term:
                    return m.group(0)
                pos = m.start()
                if _is_inside_protected(text_state[0], pos):
                    return m.group(0)
                tail = text_state[0][m.end():m.end() + 80]
                if tail.startswith(" (") or tail.startswith(" /") or tail.startswith("("):
                    return m.group(0)
                n_done[0] += 1
                return f"{m.group(0)} ({term.primary_en})"

            text_state[0] = pattern.sub(_sub, text_state[0])
            text = text_state[0]
            if n_done[0]:
                key = f"{term.ru} [{form}]" if form != term.ru else term.ru
                counts[key] = n_done[0]

    return text, counts


def enrich_file(path: Path, *, in_place: bool = False, suffix: str = ".bil.md") -> Path:
    """Обогатить файл и записать результат."""
    p = Path(path)
    src = p.read_text(encoding="utf-8")
    out, _stats = enrich_text(src)
    dst = p if in_place else p.with_suffix(suffix) if suffix.startswith(".") else p.with_name(p.name + suffix)
    dst.write_text(out, encoding="utf-8")
    return dst


# ---------------- CLI ----------------

def _main() -> int:
    """CLI: `python -m gamma.reporting.bilingual_narrator <file> [--in-place]`."""
    import argparse
    parser = argparse.ArgumentParser(
        description="F-260 Двуязычный narrative enricher (ru → ru + en). Будыка-2021 словарь."
    )
    parser.add_argument("file", help="Путь к Markdown файлу.")
    parser.add_argument("--in-place", action="store_true", help="Переписать тот же файл.")
    parser.add_argument("--suffix", default=".bil.md", help="Суффикс выходного файла.")
    parser.add_argument("--max-per-term", type=int, default=1,
                        help="Сколько раз обогащать один термин (default 1).")
    parser.add_argument("--stats", action="store_true",
                        help="Вывести статистику применённых обогащений.")
    args = parser.parse_args()

    path = Path(args.file)
    src = path.read_text(encoding="utf-8")
    glossary = load_glossary()
    out, stats = enrich_text(src, glossary, max_enrichments_per_term=args.max_per_term)

    dst = path if args.in_place else path.with_suffix(args.suffix)
    dst.write_text(out, encoding="utf-8")

    print(f"[bilingual_narrator] Output: {dst}", flush=True)
    print(f"[bilingual_narrator] Total terms enriched: {sum(stats.values())}", flush=True)
    if args.stats:
        for term, n in sorted(stats.items(), key=lambda x: -x[1]):
            print(f"  {term}: {n}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
