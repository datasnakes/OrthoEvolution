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


def test_blast_manifest_selects_exact_database_and_taxonomy() -> None:
    """Use manifest identities instead of ambiguous filename substrings."""
    client = make_client()
    metadata = [
        {
            "dbname": "refseq_rna",
            "files": [
                "ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna.00.tar.gz",
                "ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna.01.tar.gz",
            ],
        },
        {
            "dbname": "refseq_rna_index",
            "files": [
                "ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna_index.tar.gz"
            ],
        },
        {
            "dbname": "taxdb",
            "files": ["ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"],
        },
    ]

    with mock.patch.object(client, "_load_blast_metadata", return_value=metadata):
        archive_paths = client._blast_archive_paths("refseq_rna")

    assert archive_paths == [
        "/blast/db/refseq_rna.00.tar.gz",
        "/blast/db/refseq_rna.01.tar.gz",
        "/blast/db/taxdb.tar.gz",
    ]


def test_blast_manifest_rejects_unknown_database() -> None:
    """Fail before downloading when NCBI does not advertise a database."""
    client = make_client()

    with mock.patch.object(client, "_load_blast_metadata", return_value=[]):
        with pytest.raises(FileNotFoundError, match="missing"):
            client._blast_archive_paths("missing")


def test_download_does_not_replace_file_after_interruption(
    tmp_path: Path,
) -> None:
    """Keep a previous complete file when a replacement transfer fails."""
    client = make_client()
    destination = tmp_path / "database.tar.gz"
    destination.write_bytes(b"previous complete content")

    class InterruptedResponse:
        """Simulate a lost connection before a transfer completes."""

        def __enter__(self) -> "InterruptedResponse":
            return self

        def __exit__(self, *args: object) -> None:
            return None

        def read(self, _: int) -> bytes:
            raise OSError("connection lost")

    with mock.patch.object(
        client,
        "_remote_headers",
        return_value={"content-length": "100"},
    ):
        with mock.patch(
            "OrthoEvol.Tools.ftp.ncbiftp.urlopen",
            return_value=InterruptedResponse(),
        ):
            with pytest.raises(OSError, match="connection lost"):
                client._download_url("/blast/db/database.tar.gz", destination)

    assert destination.read_bytes() == b"previous complete content"
    assert not (tmp_path / "database.tar.gz.part").exists()


def test_download_rejects_checksum_mismatch(tmp_path: Path) -> None:
    """Reject a complete transfer when it does not match NCBI's checksum."""
    client = make_client()
    destination = tmp_path / "database.tar.gz"

    class BytesResponse(io.BytesIO):
        """Provide a context-managed in-memory HTTP response."""

        def __enter__(self) -> "BytesResponse":
            return self

        def __exit__(self, *args: object) -> None:
            self.close()

    with mock.patch.object(
        client,
        "_remote_headers",
        return_value={"content-length": "7"},
    ):
        with mock.patch(
            "OrthoEvol.Tools.ftp.ncbiftp.urlopen",
            return_value=BytesResponse(b"invalid"),
        ):
            with pytest.raises(OSError, match="MD5 validation failed"):
                client._download_url(
                    "/blast/db/database.tar.gz",
                    destination,
                    expected_md5="0" * 32,
                )

    assert not destination.exists()
    assert not (tmp_path / "database.tar.gz.part").exists()


def test_extract_file_rejects_parent_directory_member(tmp_path: Path) -> None:
    """Prevent an archive from writing outside the requested destination."""
    client = make_client()
    archive_path = tmp_path / "unsafe.tar.gz"
    outside_path = tmp_path.parent / "outside.txt"
    archive_member = tarfile.TarInfo("../outside.txt")
    archive_member.size = len(b"unsafe")

    with tarfile.open(archive_path, mode="w:gz") as archive:
        archive.addfile(archive_member, io.BytesIO(b"unsafe"))

    with pytest.raises(ValueError, match="Unsafe path"):
        client.extract_file(archive_path, download_path=tmp_path)

    assert not outside_path.exists()
    assert archive_path.exists()


def test_refseq_release_uses_exact_current_filename_pattern(
    tmp_path: Path,
) -> None:
    """Select all valid subparts without mixing molecule or format types."""
    client = make_client()
    release_files = [
        "vertebrate_mammalian.1.rna.gbff.gz",
        "vertebrate_mammalian.2.rna.gbff.gz",
        "vertebrate_mammalian.3.1.genomic.fna.gz",
        "vertebrate_mammalian.4.protein.gpff.gz",
    ]

    with mock.patch.object(
        client,
        "listdirectories",
        return_value=["vertebrate_mammalian"],
    ):
        with mock.patch.object(client, "listfiles", return_value=release_files):
            with mock.patch.object(client, "_read_remote_text", return_value="236"):
                with mock.patch.object(
                    client,
                    "_refseq_checksums",
                    return_value={file_name: "0" * 32 for file_name in release_files[:2]},
                ):
                    with mock.patch.object(
                        client,
                        "_download_pool",
                        return_value=[],
                    ) as download_pool:
                        client.getrefseqrelease(
                            collection_subset="vertebrate_mammalian",
                            seqtype="rna",
                            seqformat="gbff",
                            download_path=tmp_path,
                            extract=False,
                        )

    assert client.files2download == [
        "vertebrate_mammalian.1.rna.gbff.gz",
        "vertebrate_mammalian.2.rna.gbff.gz",
    ]
    download_pool.assert_called_once_with(
        [
            "/refseq/release/vertebrate_mammalian/"
            "vertebrate_mammalian.1.rna.gbff.gz",
            "/refseq/release/vertebrate_mammalian/"
            "vertebrate_mammalian.2.rna.gbff.gz",
        ],
        tmp_path,
        {
            "/refseq/release/vertebrate_mammalian/"
            "vertebrate_mammalian.1.rna.gbff.gz": "0" * 32,
            "/refseq/release/vertebrate_mammalian/"
            "vertebrate_mammalian.2.rna.gbff.gz": "0" * 32,
        },
    )


def test_refseq_release_includes_nonredundant_wp_proteins(
    tmp_path: Path,
) -> None:
    """Include NCBI's special RefSeq filename form for WP_ proteins."""
    client = make_client()
    release_files = [
        "bacteria.1.protein.gpff.gz",
        "bacteria.wp_protein.2.protein.gpff.gz",
        "bacteria.wp_protein.3.protein.faa.gz",
    ]

    with mock.patch.object(client, "listdirectories", return_value=["bacteria"]):
        with mock.patch.object(client, "listfiles", return_value=release_files):
            with mock.patch.object(client, "_read_remote_text", return_value="236"):
                with mock.patch.object(
                    client,
                    "_refseq_checksums",
                    return_value={
                        "bacteria.1.protein.gpff.gz": "0" * 32,
                        "bacteria.wp_protein.2.protein.gpff.gz": "1" * 32,
                    },
                ):
                    with mock.patch.object(
                        client,
                        "_download_pool",
                        return_value=[],
                    ):
                        client.getrefseqrelease(
                            collection_subset="bacteria",
                            seqtype="protein",
                            seqformat="gpff",
                            download_path=tmp_path,
                            extract=False,
                        )

    assert client.files2download == [
        "bacteria.1.protein.gpff.gz",
        "bacteria.wp_protein.2.protein.gpff.gz",
    ]


def test_refseq_catalog_requires_every_selected_checksum() -> None:
    client = make_client()
    catalog = "0" * 32 + "\tbacteria.1.protein.gpff.gz\n"

    with mock.patch.object(client, "_read_remote_text", return_value=catalog):
        with pytest.raises(ValueError, match="bacteria.2.protein.gpff.gz"):
            client._refseq_checksums(
                "237",
                [
                    "bacteria.1.protein.gpff.gz",
                    "bacteria.2.protein.gpff.gz",
                ],
            )
