"""Unit tests for the PHYLIP command wrapper."""

from pathlib import Path
from unittest import mock

import pexpect
import pytest

from OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip import Phylip


def build_phylip(infile: Path) -> Phylip:
    """Build only the state required by command methods."""
    phylip = object.__new__(Phylip)
    phylip.infile = infile
    phylip.phylip_log = mock.Mock()
    phylip._rename = mock.Mock()
    return phylip


def test_constructor_accepts_valid_phylip_file() -> None:
    infile = Path(__file__).parent / "test_data" / "test.phy"

    with mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip.sys.platform",
        "linux",
    ):
        phylip = Phylip(infile)

    assert phylip.infile == infile


def test_constructor_rejects_non_linux_platform() -> None:
    with mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip.sys.platform",
        "darwin",
    ):
        with pytest.raises(OSError, match="strictly for use on Linux"):
            Phylip("alignment.phy")


def test_constructor_rejects_invalid_phylip_file(tmp_path: Path) -> None:
    infile = tmp_path / "invalid.phy"
    infile.write_text("not phylip\n", encoding="utf-8")

    with mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip.sys.platform",
        "linux",
    ):
        with pytest.raises(ValueError, match="Invalid phylip format"):
            Phylip(infile)


@pytest.mark.parametrize(
    ("method_name", "arguments", "expected_outputs"),
    [
        ("dnapars", ("pars.out", "pars.tree"), ("pars.out", "pars.tree")),
        ("dnaml", ("ml.out", "ml.tree"), ("ml.out", "ml.tree")),
        ("dnadist", ("distance.out",), ("distance.out",)),
    ],
)
def test_phylip_commands_rename_outputs_and_remove_temporary_input(
    method_name: str,
    arguments: tuple[str, ...],
    expected_outputs: tuple[str, ...],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "source.phy"
    source.write_text("alignment\n", encoding="utf-8")
    phylip = build_phylip(source)
    process = mock.Mock()
    process.read.return_value = "completed"
    monkeypatch.chdir(tmp_path)

    with mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip.pexpect.spawnu",
        return_value=process,
    ):
        getattr(phylip, method_name)(*arguments)

    process.sendline.assert_called_once_with("Y\r")
    process.waitnoecho.assert_called_once_with()
    assert phylip._rename.call_args_list == [
        mock.call(source_name, destination)
        for source_name, destination in zip(
            ("outfile", "outtree"),
            expected_outputs,
            strict=False,
        )
    ]
    assert not (tmp_path / "infile").exists()


@pytest.mark.parametrize(
    ("method_name", "arguments"),
    [
        ("dnapars", ("pars.out", "pars.tree")),
        ("dnaml", ("ml.out", "ml.tree")),
        ("dnadist", ("distance.out",)),
    ],
)
def test_phylip_commands_log_eof_and_remove_temporary_input(
    method_name: str,
    arguments: tuple[str, ...],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "source.phy"
    source.write_text("alignment\n", encoding="utf-8")
    phylip = build_phylip(source)
    process = mock.Mock()
    process.waitnoecho.side_effect = pexpect.EOF("finished")
    process.read.return_value = "diagnostic"
    monkeypatch.chdir(tmp_path)

    with mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.Phylip.phylip.pexpect.spawnu",
        return_value=process,
    ):
        getattr(phylip, method_name)(*arguments)

    phylip.phylip_log.error.assert_called_once_with("diagnostic")
    phylip._rename.assert_not_called()
    assert not (tmp_path / "infile").exists()
