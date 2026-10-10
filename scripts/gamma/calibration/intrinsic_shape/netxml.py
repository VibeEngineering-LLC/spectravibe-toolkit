"""Разностный спектр «образец минус фон» после регистрации по форме -> два XML для BecqMoni.
Перенесено из прототипа labr3_selfcal/make_net_xml.py:123-143 (write_xml), 2026-10-10: текст исходного XML образца
сохраняется, заменяются имя, примечание, шкала, Valid/TotalPulseCount, DataPoint; добавляется блок SqrtFwhmCalibration.
DataPoint -- только целые (в BecqMoni Spectrum -- int[]).
"""
import re
import numpy as np

from .fwhm import fwhm_block, insert_block


def write_xml(src, path, data, signed, C, H, k, fwc, fw_peaks, name, date, times=None):
    """src -- XML образца (шаблон текста); data -- целые счёты на шкале фона; C -- шкала E(ch) фона;
    H -- отображение h; fwc -- коэффициенты ПШПВ; fw_peaks -- [(канал, ПШПВ)]; signed -- файл «со знаком»;
    times -- None (времена образца) или (LiveTime, MeasurementTime) выходного файла (правка Claude, раунд 4)."""
    data = np.asarray(np.round(data), dtype=np.int64)   # DataPoint только целые
    t = open(src, encoding="utf-8", newline="").read()
    if times is not None:   # первое вхождение -- спектр образца (блок фона ниже удаляется)
        t = re.sub(r"<LiveTime>[^<]*</LiveTime>", "<LiveTime>%r</LiveTime>" % float(times[0]), t, count=1)
        t = re.sub(r"<MeasurementTime>[^<]*</MeasurementTime>", "<MeasurementTime>%d</MeasurementTime>" % int(round(times[1])), t, count=1)
    nl = "\r\n" if "\r\n" in t else "\n"
    t = re.sub(r"(<SampleInfo>\s*<Name>)[^<]*(</Name>)", lambda m: m.group(1) + name + m.group(2), t, count=1)
    note = "Отображение h(x)=%.7g+%.7g*x%+.4g*x^2 (по форме), k=%.6f, шкала по форме, %s" % (H[0], H[1], H[2], k, date)
    t = re.sub(r"<Note />|<Note>[^<]*</Note>", lambda m: "<Note>%s</Note>" % note, t, count=1)
    # правка Claude (раунд 4, п.7): разность уже вычтена -- фон из файла убирается целиком (имя и блок спектра),
    # иначе BecqMoni в режиме вычитания фона вычтет его второй раз. Без них BackgroundEnergySpectrum = null
    # (значение по умолчанию поля, ResultData.cs:429), имя фона = "" (ResultData.cs:423); null проверяется везде.
    t = re.sub(r"[ \t]*<BackgroundSpectrumFile\s*/>\r?\n|[ \t]*<BackgroundSpectrumFile>[^<]*</BackgroundSpectrumFile>\r?\n", "", t, count=1)
    t = re.sub(r"[ \t]*<BackgroundEnergySpectrum>.*?</BackgroundEnergySpectrum>\r?\n|[ \t]*<BackgroundEnergySpectrum\s*/>\r?\n",
               "", t, count=1, flags=re.S)
    t = re.sub(r"<FormatVersion>120920</FormatVersion>", "<FormatVersion />", t, count=1)
    t = re.sub(r"<PolynomialOrder>\d+</PolynomialOrder>", "<PolynomialOrder>%d</PolynomialOrder>" % (len(C) - 1), t, count=1)
    t = re.sub(r"(<Coefficients>).*?(</Coefficients>)", lambda m: m.group(1) + nl + "".join(
        "            <Coefficient>%.17g</Coefficient>%s" % (c, nl) for c in C) + "          " + m.group(2), t, count=1, flags=re.S)
    # Valid/Total: обнулённый -- сумма содержимого (как в прототипе); «со знаком» -- Valid = сумма содержимого
    # (вторая проверка: было сумма положительных), Total = Valid + (Total - Valid исходника), если Total > Valid.
    valid = int(data.sum())
    sv = re.search(r"<ValidPulseCount>(\d+)</ValidPulseCount>", t)
    st = re.search(r"<TotalPulseCount>(\d+)</TotalPulseCount>", t)
    extra = max(0, int(st.group(1)) - int(sv.group(1))) if (signed and sv and st) else 0
    t = re.sub(r"<ValidPulseCount>-?\d+</ValidPulseCount>", "<ValidPulseCount>%d</ValidPulseCount>" % valid, t, count=1)
    t = re.sub(r"<TotalPulseCount>-?\d+</TotalPulseCount>", "<TotalPulseCount>%d</TotalPulseCount>" % (valid + extra), t, count=1)
    pts_xml = nl.join("          <DataPoint>%d</DataPoint>" % v for v in data)
    t = insert_block(t, fwhm_block(nl, fwc, fw_peaks, C))
    t = re.sub(r"(<Spectrum>)\r?\n.*?(\r?\n\s*</Spectrum>)", lambda m: m.group(1) + nl + pts_xml + m.group(2), t, count=1, flags=re.S)
    open(path, "w", encoding="utf-8", newline="").write(t)
