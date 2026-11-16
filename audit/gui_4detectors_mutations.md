# Мутационная поверка тестов окна: 4 детектора

Тесты: `apps/intrinsic-shape-gui/tests/test_gui_detectors.py`; прогон: `python apps/intrinsic-shape-gui/tests/gui_mutations.py`.

Исправная версия: 11 passed in 4.77s

| Мутация | Что портит | Файл | Итог | Покраснели тесты |
|---|---|---|---|---|
| M1_profiles | список профилей: только La-138 | backend_nai.py | убита | `test_about_and_version`, `test_autodetect_switches_to_la_profile`, `test_broken_xml_friendly_error`, `test_cancel_nai_run`, `test_combo_filled_from_profile_layer`, `test_foreign_method_tabs_hidden`, `test_mismatch_is_warning_not_refusal`, `test_nai_hint_without_match`, `test_nai_run_writes_experimental_xml`, `test_no_path_to_foreign_command`, `test_selection_persisted` |
| M2_tabs | скрытие вкладок чужого метода | main_window.py | убита | `test_foreign_method_tabs_hidden` |
| M3_banner | пометка «экспериментально» на вкладке | nai_panel.py | убита | `test_nai_run_writes_experimental_xml` |
| M3b_summary | пометка «экспериментально» в итоге | backend_nai.py | убита | `test_nai_run_writes_experimental_xml` |
| M4_guard | выбор метода: проверка метода команды | main_window.py | убита | `test_no_path_to_foreign_command` |
| M4b_naimethod | выбор метода: запуск NaI как La-138 | main_window.py | убита | `test_cancel_nai_run`, `test_nai_run_writes_experimental_xml`, `test_no_path_to_foreign_command` |
| M5_detect | автоопределение: прибор не распознаётся | backend_nai.py | убита | `test_autodetect_switches_to_la_profile`, `test_mismatch_is_warning_not_refusal` |
| M5b_autoswitch | автоопределение: профиль не переключается | main_window.py | убита | `test_autodetect_switches_to_la_profile` |
| M6_settings | сохранение выбора профиля | main_window.py | убита | `test_selection_persisted` |
| M6b_autoflag | сохранение флага автовыбора | main_window.py | убита | `test_selection_persisted` |
| M7_cancel | отмена | backend_nai.py | убита | `test_cancel_nai_run` |
| M8_mismatch | предупреждение о несовпадении (распознан другой) | backend_nai.py | убита | `test_mismatch_is_warning_not_refusal` |
| M8b_mismatch_win | предупреждение о несовпадении в окне | main_window.py | убита | `test_mismatch_is_warning_not_refusal` |
| M9_peaks_nuclides | режим по линиям: явный список нуклидов | backend_nai.py | убита | `test_nai_run_writes_experimental_xml` |

Итого: мутаций 14, убито 14, выжило 0.
