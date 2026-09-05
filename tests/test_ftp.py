"""Tests for FTP behaviors that do not require a network connection."""

import io
import tarfile
from pathlib import Path
from unittest import mock

import pytest

from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Tools.ftp.baseftp import BaseFTPClient
from OrthoEvol.Tools.ftp.ncbiftp import NcbiFTPClient


def make_client() -> NcbiFTPClient:
    """Build an NCBI client without opening a network connection."""
    client = object.__new__(NcbiFTPClient)
    client.email = "test@example.org"
    client.cpus = 2
    client._timeout = 30.0
    client.files2download = []
    client.refseqrelease_path = "/refseq/release/"
    client.refseq_release_number_path = "/refseq/release/RELEASE_NUMBER"
    client.ncbiftp_log = mock.Mock()
    return client


def test_windowmasker_download_remains_explicitly_unsupported() -> None:
    """Preserve the deprecated API boundary without opening an FTP connection."""
    client = make_client()

    with pytest.raises(
        OrthoEvolDeprecationWarning,
        match="WindowMasker downloads are no longer supported",
    ):
        client.getwindowmaskerfiles([9606], "/tmp")


def test_keepalive_transfer_preserves_same_named_local_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Discard keepalive bytes without touching a caller-owned local file."""
    client = object.__new__(BaseFTPClient)
    client.ftp = mock.Mock()
    client.ftp.pwd.return_value = "/remote/current"
    client.ftp.retrbinary.side_effect = (
        lambda _command, callback: callback(b"remote content")
    )
    local_file = tmp_path / "README.ftp"
    local_file.write_text("local content", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    client._filetransfer("README.ftp")

    assert local_file.read_text(encoding="utf-8") == "local content"
    assert client.ftp.cwd.call_args_list == [
        mock.call("/"),
        mock.call("/remote/current"),
    ]
