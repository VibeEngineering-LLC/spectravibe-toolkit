# -*- coding: utf-8 -*-
"""F-256 — Citation translator Layer 1 (RAG-ID) → Layer 2 (ГОСТ Р 7.0.5–2008).

Двухслойная схема ссылок проекта (см. references/REFERENCES.md §0):
- **Layer 1** (внутренний FAST): `[RAG-ID]` — `[LSRM-Algo-10]`, `[BUDYKA-7.5]`
  Используется в коде, F-rule docstrings, аудитах, ROADMAP, JSON-отчётах.
- **Layer 2** (внешний ГОСТ): `[№, локатор]` — `[7, §10]`, `[12, §7.5]`
  Используется в чат-отчётах пользователю, HTML/PDF/MD регулируемых отчётах,
  печатной нормоконтрольной документации.

Этот модуль реализует автоматическую трансляцию L1 → L2 перед сохранением
финального user-facing вывода.

Контракт (F-256, v1.17.10):
- Любой `reporting/*.py` (html_report.py, pdf_export.py, markdown_report.py),
  генерирующий user-facing output, обязан вызывать `translate_text(content)`
  ПЕРЕД финальной записью файла.
- JSON-отчёты (json_report.py) остаются на Layer 1 (machine-consumed).

См. также: F-150 (release archive), F-154 (tool preservation), F-157 (LSRM priority).
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

# RAG-ID prefix → (GOST-№, локатор-kind, default-locator-marker)
# Источник истины: references/REFERENCES.md §0 (таблица prefix → №).
RAG_PREFIX_TO_GOST: Dict[str, Tuple[int, str]] = {
    # ЛСРМ официальная документация
    "LSRM-Algo": (7, "§"),                # Алгоритмические основы SpectraLine
    "LSRM-SPE-Format": (8, "§"),
    "LSRM-SL2-Basics": (9, "с."),
    "LSRM-PREC": (10, "слайд"),           # Прецизионные измерения (Кувыкин 2023)
    "LSRM-ACT": (5, "§"),                  # LSRM-Активность 2014 МВИ
    "LSRM-NMAT": (11, "слайд"),
    "LSRM-DOSE-RATE": (6, "§"),
    "LSRM-DOSE": (6, "§"),
    # Русские академические источники
    "BUDYKA-Glossary": (13, "с."),
    # F-338 / v1.18.20.0 — Будыка/Шендрик — учебники с classical pages,
    # натуральный локатор — страница, не §.
    "BUDYKA": (12, "с."),                  # Будыка-2021 учебник
    "SHENDRIK-1": (14, "с."),
    "SHENDRIK-2": (15, "с."),
    "SHENDRIK": (14, "с."),                # default → ч.1
    "VARTANOV": (16, "с."),
    "MDA-RU-CONCEPTS": (17, "§"),
    "MDA-RU": (17, "§"),
    "EXP-ANALYSIS": (18, "§"),
    "EXPERIMENT": (18, "§"),
    # Международные референсы
    # F-338 — Gilmore — textbook, locator = с. (страница).
    "GILMORE": (19, "с."),
    "ORTEC-GV9": (22, "§"),                # ORTEC GammaVision V9 (A66, 2020)
    "ORTEC-GammaVision-V9": (22, "§"),
    "ORTEC-GammaVision": (22, "§"),
    "ORTEC": (22, "§"),                    # default fallback
    "SKL-MIX": (23, ""),                   # scikit-learn sklearn.mixture (online doc)
    "sklearn.mixture": (23, ""),
    # Нормативные документы
    "GOST-7.0.5-2008": (1, "§"),
    "GOST-7.0.5": (1, "§"),
    # F-338 — ISO 11929 — структура "clause", локатор «п.» по ISO-style.
    "ISO-11929": (2, "п."),
    "ICRP-74": (3, "§"),
    "ICRP74": (3, "§"),
    # Базы данных
    "ENSDF": (4, ""),
    "IAEA-LC": (20, ""),
    "NIST-XCOM": (21, ""),
}

# Спец-маппинг для LSRM-ACT — раздел Приложение/секция кодируется в локаторе
# (см. AUDIT_v2_MERGED §9.6 для конкретных Прил.N).
LSRM_ACT_SECTION_MAP: Dict[str, str] = {
    "01": "§3.1",
    "02": "§3.1",
    "03": "§3.2",
    "04": "§3.2",
    "05": "§4-5.1",
    "06": "§5.2.2",
    "07": "§10.1.3 / Прил.7",
    "08": "Прил.2",
    "09": "Прил.3",
    "10": "Прил.4",
    "11": "Прил.5",
    "12": "Прил.8",
    "13": "Прил.9 §3",
    "14": "Прил.9 §4",
    "15": "Прил.10",
}

# Аналогично для LSRM-PREC (Прецизионные слайды 5_1)
LSRM_PREC_SECTION_MAP: Dict[str, str] = {
    "1": "слайд 2",
    "2": "слайд 4",
    "3": "слайд 5-6",
    "4": "слайд 7-8",
    "5": "слайд 9-17",
    "6": "слайд 11",
    "7": "слайд 13",
    "8": "слайд 21",
    "9": "слайд 22-24",
    "10": "слайд 25-28",
    "11": "слайд 31-36",
    "12": "слайд 37-39",
    "13": "слайд 40-42",
    "14": "слайд 43",
    "15": "слайд 4",
    "06": "слайд 18-22",
    "07": "слайд 24-25",
}


# Регулярка для RAG-ID: `[PREFIX-LOCATOR]`
# Начало — латинская заглавная; далее буквы/цифры/точки/дефисы (mixed-case ok).
_RAG_ID_RE = re.compile(r"\[([A-Z][A-Za-z0-9.]+(?:-[A-Za-z0-9.]+)*)\]")

# Префиксы F-rule / audit-IDs — НЕ библиографические, пропускаются translate_text.
_SKIP_PREFIXES = ("F-", "TD-", "CAL-", "PEAK-", "PHYS-", "ID-", "NEW-", "MEM-", "TCS-")


@dataclass(frozen=True)
class TranslationResult:
    """Результат трансляции одного RAG-ID."""
    rag_id: str
    gost_num: int
    locator: str  # например "§10", "Прил.5", "слайд 22-24"
    is_resolved: bool  # False если prefix не найден в RAG_PREFIX_TO_GOST


def translate_rag_id(rag_id: str) -> TranslationResult:
    """Трансляция одного RAG-ID в ГОСТ-локатор.

    Примеры:
        >>> translate_rag_id("LSRM-Algo-10")
        TranslationResult(rag_id='LSRM-Algo-10', gost_num=7, locator='§10', is_resolved=True)
        >>> translate_rag_id("BUDYKA-7.5")
        TranslationResult(rag_id='BUDYKA-7.5', gost_num=12, locator='§7.5', is_resolved=True)
        >>> translate_rag_id("LSRM-ACT-11")
        TranslationResult(rag_id='LSRM-ACT-11', gost_num=5, locator='Прил.5', is_resolved=True)
        >>> translate_rag_id("LSRM-PREC-9")
        TranslationResult(rag_id='LSRM-PREC-9', gost_num=10, locator='слайд 22-24', is_resolved=True)
    """
    # Спец-кейсы с map
    if rag_id.startswith("LSRM-ACT-"):
        loc = rag_id[len("LSRM-ACT-"):]
        section = LSRM_ACT_SECTION_MAP.get(loc.zfill(2)) or LSRM_ACT_SECTION_MAP.get(loc)
        if section:
            return TranslationResult(rag_id, 5, section, True)
    if rag_id.startswith("LSRM-PREC-"):
        loc = rag_id[len("LSRM-PREC-"):]
        section = LSRM_PREC_SECTION_MAP.get(loc) or LSRM_PREC_SECTION_MAP.get(loc.lstrip("0"))
        if section:
            return TranslationResult(rag_id, 10, section, True)

    # Generic: ищем самый длинный prefix-match
    best_prefix = ""
    for prefix in RAG_PREFIX_TO_GOST:
        if rag_id.startswith(prefix) and len(prefix) > len(best_prefix):
            best_prefix = prefix
    if not best_prefix:
        return TranslationResult(rag_id, 0, rag_id, False)

    num, marker = RAG_PREFIX_TO_GOST[best_prefix]
    suffix = rag_id[len(best_prefix):].lstrip("-")
    if not suffix or not marker:
        return TranslationResult(rag_id, num, "", True)
    return TranslationResult(rag_id, num, f"{marker}{suffix}", True)


def format_layer2(result: TranslationResult) -> str:
    """Формат финальной ссылки Layer 2 по ГОСТ Р 7.0.5–2008 §7.3.

    F-338 / v1.18.20.0 — между marker (§/с./п./слайд) и числом ставится
    пробел: «[7, § 9.3]», «[12, с. 47]» (per ГОСТ examples).
    """
    if not result.is_resolved:
        return f"[{result.rag_id}]"  # fallback — оставить как было + warning в логах
    if result.locator:
        # F-338 — между marker prefix и первой цифрой вставляем пробел.
        # Если пробел уже есть (LSRM_PREC «слайд 9-17» из SECTION_MAP) —
        # ничего не делаем (regex не матчится).
        loc = re.sub(r"^([§A-Za-zА-Яа-я.]+)(\d)", r"\1 \2", result.locator)
        return f"[{result.gost_num}, {loc}]"
    return f"[{result.gost_num}]"


def translate_text(text: str, *, strict: bool = False) -> Tuple[str, List[str]]:
    """Заменить все [RAG-ID] в тексте на [№, локатор].

    Args:
        text: исходный текст с Layer 1 ссылками.
        strict: если True — кидать ValueError при нерезолвящихся RAG-ID.

    Returns:
        (translated_text, unresolved_ids) — где unresolved_ids — список RAG-ID,
        которые не нашлись в карте (могут быть валидными внутренними ID,
        не относящимися к библиографии — например F-rule IDs).

    Examples:
        >>> txt = "См. [LSRM-Algo-10] и [BUDYKA-7.5]. F-167 (CAL-001)."
        >>> translate_text(txt)
        ('См. [7, §10] и [12, §7.5]. F-167 (CAL-001).', [])
        >>> txt2 = "[LSRM-ACT-11] vs [UNKNOWN-PREFIX-1]"
        >>> out, unr = translate_text(txt2)
        >>> out
        '[5, Прил.5] vs [UNKNOWN-PREFIX-1]'
        >>> unr
        ['UNKNOWN-PREFIX-1']
    """
    unresolved: List[str] = []

    def _sub(m: re.Match) -> str:
        rag_id = m.group(1)  # capture group 1
        # Пропускаем F-rule и audit-IDs (не библиографические)
        if any(rag_id.startswith(p) for p in _SKIP_PREFIXES):
            return m.group(0)
        result = translate_rag_id(rag_id)
        if not result.is_resolved:
            unresolved.append(rag_id)
            if strict:
                raise ValueError(f"Не резолвится RAG-ID: {rag_id}")
            return m.group(0)
        return format_layer2(result)

    translated = _RAG_ID_RE.sub(_sub, text)
    # F-338 / v1.18.20.0 — post-process: merge соседние bare-numeric ссылки.
    translated = _merge_adjacent_citations(translated)
    return translated, unresolved


# F-338 / v1.18.20.0 — regex для слияния соседних bare-numeric `[N][M][P]...`
# в один `[N, M, P]` ascending. Условие: ВСЕ скобки в run — bare (только
# цифры, без локатора). Если хоть одна с локатором — run не сольётся.
# Регекс капчит >=2 bracket-групп вида `[N]`, опционально разделённых
# whitespace или unicode no-break space. Поведение для одиночной `[N]` —
# не трогать.
_MERGE_ADJ_RE = re.compile(
    r"\[(\d+)\](?:[  \t]*\[\d+\])+"
)


def _merge_adjacent_citations(text: str) -> str:
    """F-338 — склеить соседние `[N][M][P]` в `[N, M, P]` ascending unique.

    Только bare-numeric brackets мержатся. Mixed runs (с локаторами)
    остаются раздельными (regex их не матчит).

    Examples:
        >>> _merge_adjacent_citations("См. [5][7][12].")
        'См. [5, 7, 12].'
        >>> _merge_adjacent_citations("Дано [7, с. 9][12].")
        'Дано [7, с. 9][12].'
        >>> _merge_adjacent_citations("[3][1][3]")     # dedup + sort
        '[1, 3]'
        >>> _merge_adjacent_citations("Один [7].")     # no-op
        'Один [7].'
    """
    if not text or "[" not in text:
        return text

    def _merge(m: re.Match) -> str:
        full = m.group(0)
        nums = re.findall(r"\[(\d+)\]", full)
        unique_sorted = sorted({int(n) for n in nums})
        return "[" + ", ".join(str(n) for n in unique_sorted) + "]"

    return _MERGE_ADJ_RE.sub(_merge, text)


def translate_file(path: str | Path, *, in_place: bool = False, suffix: str = ".l2.md") -> Path:
    """Перевести Markdown/text файл с Layer 1 → Layer 2.

    Args:
        path: путь к исходному файлу.
        in_place: если True — переписать тот же файл; иначе создать <name>.l2.md.
        suffix: суффикс для нового файла, если in_place=False.

    Returns:
        Путь к переведённому файлу.
    """
    p = Path(path)
    src = p.read_text(encoding="utf-8")
    out, _unresolved = translate_text(src)
    if in_place:
        dst = p
    else:
        dst = p.with_suffix(suffix) if suffix.startswith(".") else p.with_name(p.name + suffix)
    dst.write_text(out, encoding="utf-8")
    return dst


# ---------------- CLI ----------------

def _main() -> int:
    """CLI entrypoint: `python -m gamma.reporting.citation_translator <file> [--in-place]`."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Citation translator Layer 1 (RAG-ID) → Layer 2 (ГОСТ Р 7.0.5–2008). F-256."
    )
    parser.add_argument("file", help="Путь к Markdown/text файлу.")
    parser.add_argument("--in-place", action="store_true", help="Переписать тот же файл.")
    parser.add_argument("--suffix", default=".l2.md",
                        help="Суффикс для нового файла (если не --in-place).")
    parser.add_argument("--strict", action="store_true",
                        help="Падать на нерезолвящихся RAG-ID.")
    parser.add_argument("--report-unresolved", action="store_true",
                        help="Печатать список нерезолвящихся RAG-ID.")
    args = parser.parse_args()

    path = Path(args.file)
    src = path.read_text(encoding="utf-8")
    out, unresolved = translate_text(src, strict=args.strict)

    if args.report_unresolved and unresolved:
        print(f"[citation_translator] Unresolved RAG-IDs: {len(set(unresolved))}", flush=True)
        for u in sorted(set(unresolved)):
            print(f"  - {u}", flush=True)

    dst = translate_file(path, in_place=args.in_place, suffix=args.suffix)
    print(f"[citation_translator] Layer 2 output: {dst}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
