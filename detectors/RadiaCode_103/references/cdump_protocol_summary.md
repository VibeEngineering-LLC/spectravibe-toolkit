# cdump radiacode — короткий operator-facing конспект

**Дата фиксации**: 2026-06-14
**Источник**: https://github.com/cdump/radiacode (master, 2026-06-14)
**Полная техническая справка**: `<USERPROFILE>\.claude\skills\radiacode-ble\references\cdump_radiacode_python.md` и `protocol.md`. Здесь — operator-side выжимка ≈ 1-2 KB про практический workflow.
**Связка с нашим pipeline**: RC-103 XML (`Spectrum *.xml` из RadiaCode Studio) парсится `gamma.io.becqmoni_xml`/`atomspectra_xml` напрямую → можно прогонять через `scripts/run_plan_a.py` без cdump.

## 1. Зачем cdump оператору этого проекта

`cdump/radiacode` — python-библиотека для **прямого** общения с RC-101/102/103/110 через USB или BLE (минуя RadiaCode Studio). Полезна для трёх задач:

1. **Выгрузка спектра без UI** — `examples/basic.py` сохраняет XML, совместимый с нашим pipeline'ом (BecqMoni `ResultDataFile`).
2. **Real-time мониторинг** — stream cps/μSv/h в собственный скрипт (для autonomic monitor'а, telegram-бота, MQTT, push в Home Assistant).
3. **Reset accumulation + установка времени** — без захода в родное приложение, чтобы начать чистый длинный набор фона/калибровки.

## 2. Минимальный пример выгрузки спектра в XML

```python
from radiacode import RadiaCode
from radiacode.transports.usb import Usb   # или Bluetooth

rc = RadiaCode(transport=Usb())               # требует USB-кабель + driver ST (VID 0x0483/PID 0xF123)
rc.set_local_time()                           # синк часов прибора
print(rc.fw_version())                        # должно быть ≥ 4.8.0
print(rc.serial_number())                     # ASCII RC-103-000NNN

spec = rc.spectrum()                          # тек. спектр (с накопления)
print(f"duration={spec.duration}s, counts={sum(spec.counts)}")

# Сохранение в формате, который читает наш pipeline:
import xml.etree.ElementTree as ET
# ... см. examples/basic.py в cdump для готового шаблона ...
```

**FW ≥ 4.8 обязателен** для cdump master. Если на устройстве FW < 4.8 — обновить через родное приложение RadiaCode, либо использовать `radiacode==0.2.2` (legacy без проверки).

## 3. Single-central rule

RC-устройства поддерживают **только одно** активное BLE-соединение. Перед подключением через cdump:

1. На iPhone/Android: отключить RadiaCode приложение **или** Bluetooth целиком.
2. На Windows: убедиться что RadiaCode Studio не открыта (она съедает USB-CDC порт).
3. Если соединение зависло — `power-cycle` устройство (выкл/вкл кнопкой).

При нарушении single-central rule cdump падает с `BleakError: Failed to connect` или (USB) `serial.SerialException`.

## 4. Reset accumulation — чистый длинный набор

```python
rc.spectrum_reset()   # обнулить накопление спектра (foreground)
# rc.dose_reset()     # обнулить дозу (если нужно)
# теперь набирать сколько надо, потом rc.spectrum() читать
```

Применение для нашего pipeline:
- Чистый фон Gamma-1C → reset → набирать N часов → выгрузить XML → `run_plan_a.py GAMMA_BG=...`.
- Чистый sample → reset → измерить → XML → `run_plan_a.py GAMMA_SAMPLE=... GAMMA_BG=...`.

## 5. Real-time data stream

```python
for data in rc.data_buf():
    # data — поток объектов RealTimeData / DoseRateDB / RareData / etc.
    # См. cdump_radiacode_python.md секция 5 для полной таблицы.
    if hasattr(data, 'dose_rate'):
        print(f"{data.dose_rate} µSv/h")
```

Для нашего pipeline эти данные **не используются** напрямую (мы анализируем накопленный спектр), но могут быть полезны для:
- Live-monitor виджета на дашборде.
- Отслеживания дозы оператора во время длинной поверки.
- Cross-валидации: суммарный dose из data_buf за период должен ≈ соответствовать ε(E)·Σ(counts·E) спектра.

## 6. Установка времени прибора

```python
rc.set_local_time()                           # быстро, по локальной TZ
# или
from datetime import datetime
rc.set_time(datetime(2026, 6, 14, 12, 0, 0))  # явное UTC время
```

**Зачем это в нашем pipeline**: `<StartTime>` / `<EndTime>` в выгруженном XML формируются по часам прибора. Если прибор отстаёт — даты в нашем reports будут неверные. Перед длинным набором (поверка, фон, sample) — `set_local_time()` обязателен.

## 7. Cross-reference с проектом

- **Парсер RC-103 XML**: `scripts/gamma/io/atomspectra_xml.py` (reader), `becqmoni_xml.py` (writer). Узнаёт по корневому тегу `<ResultDataFile>` через `format_registry`.
- **Operator workflow без cdump**: RadiaCode Studio → File → Export → XML → передать через `GAMMA_SAMPLE=...` в `run_plan_a.py`. Подтверждено 2026-06-13 на `Spectrum 13-06-2026.xml` (5.3 сек end-to-end).
- **Operator workflow с cdump** (опционально): `examples/basic.py` cdump'а + наш `run_plan_a.py` → можно полностью обойти RadiaCode Studio в автоматизированном pipeline.

## 8. Известные ограничения cdump

- **FW threshold 4.8**: ниже — не работает (либо legacy 0.2.2).
- **BLE chunk = 18 байт hard-coded**: длинные ответы прибор разбивает на чанки, cdump склеивает по `<i> length` префиксу.
- **RD_VIRT_STRING НDA новой FW** — добавлен trailing `0x00` маркер. Cdump обрабатывает через HACK в `_read_data` (см. `protocol.md` секция 11.7).
- **AlarmLimits единицы** — count_rate в counts/10s, dose в µR (не µSv/h!). При записи через cdump нужно конвертировать.

## 9. Безопасность и приватность

Серийник прибора (`RC-103-000NNN`) — формально приватный (можно сопоставить с покупкой), но **не критичный**. MAC (`AddressBLE`) — тоже идентификатор. Для публичных скиллов / репозиториев / отчётов — замены на `RC-103-XXXXXX` и `XX:XX:XX:XX:XX:XX`. В operator-internal reports — оставляем как есть.

## 10. Источники для дальнейшего чтения

- **Полный wire-level протокол**: `<USERPROFILE>\.claude\skills\radiacode-ble\references\protocol.md` (комбинированный — оригинальный + cdump refresh секции 11).
- **API класса RadiaCode**: `<USERPROFILE>\.claude\skills\radiacode-ble\references\cdump_radiacode_python.md` (методы, типы, dataclasses).
- **Карта моделей RC-1xx**: `<USERPROFILE>\.claude\skills\radiacode-ble\references\known_devices.md` (serial-prefix, capability matrix).
- **Reverse-engineering методология**: skill `device-protocol-re` (BLE/USB sniffing, decoding).

cdump-репо — основной reference. Здесь — выжимка для оператора этого проекта.