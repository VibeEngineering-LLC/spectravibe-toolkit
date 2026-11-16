"""Общий слой профилей: NaI-профили (метод nai, экспериментально) рядом с LaBr; команды чужого метода -- код 2."""
import pytest
from gamma.calibration.intrinsic_shape import profiles, cli

ALL = {"labr_atomnano5pro", "labr_as2pro", "nai_gamma1s", "nai_gs2020"}


def test_all_methods_ids_and_aliases():
    assert set(profiles.ids(all_methods=True)) == ALL
    assert set(profiles.ids()) == {"labr_atomnano5pro", "labr_as2pro"}
    assert profiles.get("gamma1s")["id"] == "nai_gamma1s" and profiles.get("gs2020")["id"] == "nai_gs2020"
    for pid in ("nai_gamma1s", "nai_gs2020"):
        p = profiles.get(pid)
        assert p["method"] == "nai" and p["status"] == "experimental" and "экспериментально" in p["title"]


@pytest.mark.parametrize("argv", [["register", "--detector", "nai_gamma1s"], ["nai-shape", "--detector", "labr_as2pro", "--out-dir", "."]])
def test_cli_rejects_foreign_method(argv, tmp_path, capsys):
    argv = argv[:1] + ["--bg", str(tmp_path / "b.xml"), "--sample", str(tmp_path / "s.xml")] + argv[1:]
    assert cli.main(argv) == 2
    assert "использует метод" in capsys.readouterr().err
    profiles.activate(None)
