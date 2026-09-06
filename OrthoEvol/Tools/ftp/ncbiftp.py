"""Download current NCBI BLAST and RefSeq data."""

import gzip
import hashlib
import json
import os
import re
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from email.utils import parsedate_to_datetime
from ftplib import error_perm
from pathlib import Path, PurePosixPath
from time import perf_counter
from typing import Any, NoReturn, Sequence
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from tqdm import tqdm

from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Tools.ftp.baseftp import BaseFTPClient
from OrthoEvol.Tools.logit import LogIt


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI directories over FTP and download data over HTTPS."""

    _ncbi_host = "ftp.ncbi.nlm.nih.gov"
    _blast_metadata_path = "/blast/db/blastdb-metadata-1-1.json"
    _transfer_chunk_size = 1024 * 1024

    def __init__(
        self,
        email: str,
        max_workers: int = 4,
        **kwargs: Any,
    ) -> None:
        """Initialize an anonymous NCBI connection.

        :param email: Contact email used as the anonymous FTP password.
        :param max_workers: Maximum number of concurrent HTTPS transfers.
        :param kwargs: Additional arguments for :class:`BaseFTPClient`.
        """
        if max_workers < 1:
            raise ValueError("max_workers must be at least 1.")

        self.email = email
        self.cpus = max_workers
        super().__init__(
            self._ncbi_host,
            user="anonymous",
            password=email,
            **kwargs,
        )

        self._datefmt = "%m-%d-%Y@%I:%M:%S-%p"
        self._date = datetime.now().strftime(self._datefmt)
        self.blastpath = "/blast/"
        self.blastdb_path = "/blast/db/"
        self.blastdbv4_path = "/blast/db/v4/"
        # NCBI serves version 5 databases from the main BLAST DB directory.
        self.blastdbv5_path = self.blastdb_path
        self.blastfasta_path = "/blast/db/FASTA/"
        self.refseqrelease_path = "/refseq/release/"
        self.refseq_release_number_path = "/refseq/release/RELEASE_NUMBER"
        self.windowmasker_path = "/blast/windowmasker_files/"

        # Retain these public attributes for callers that inspect selections.
        self._taxdb = ["taxdb.tar.gz", "taxdb.tar.gz.md5"]
        self.blastdbs: list[str] = []
        self.blastfastadbs: list[str] = []
        self.files2download: list[str] = []
        self.refseqreleasedbs: list[str] = []
        self.refseqrelease_seqtypes: list[str] = []
        self.refseqrelease_filetypes: list[str] = []

        self.ncbiftp_log = LogIt().default(logname="NCBI-FTP", logfile=None)

    @classmethod
    def _pathformat(cls, path: str) -> None:
        """Validate an absolute FTP directory path.

        :param path: FTP path ending in a slash.
        """
        if not path.startswith("/") or not path.endswith("/"):
            raise ValueError("FTP paths must start and end with '/'.")

    @classmethod
    def _archive(
        cls,
        archive_name: str,
        folder2archive: str | Path,
        archive_type: str = "gztar",
    ) -> Path:
        """Archive a directory without changing the process working directory.

        :param archive_name: Output archive base name.
        :param folder2archive: Directory to archive.
        :param archive_type: Format accepted by :func:`shutil.make_archive`.
        :return: Path to the created archive.
        """
        source = Path(folder2archive)
        if not source.is_dir():
            raise NotADirectoryError(source)

        archive_base = source.parent / archive_name
        archive_path = shutil.make_archive(
            str(archive_base),
            archive_type,
            root_dir=source.parent,
            base_dir=source.name,
        )
        return Path(archive_path)

    def walk(self, path: str) -> tuple[list[str], list[str]]:
        """List immediate directories and files in an NCBI FTP path.

        :param path: Absolute FTP directory path.
        :return: Directory names followed by file names.
        """
        self._pathformat(path)
        try:
            self.ftp.cwd(path)
        except error_perm as error:
            self.ncbiftp_log.info(
                "Unable to access %s from %s: %s",
                path,
                self.ftp.pwd(),
                error,
            )
            return [], []

        directories: list[str] = []
        files: list[str] = []
        try:
            for name, facts in self.ftp.mlsd():
                entry_type = facts.get("type")
                if entry_type == "dir":
                    directories.append(name)
                elif entry_type == "file":
                    files.append(name)
        except (AttributeError, error_perm):
            # LIST is retained for older FTP servers that do not implement MLSD.
            listing: list[str] = []
            self.ftp.retrlines("LIST", listing.append)
            for row in listing:
                fields = row.split()
                if not fields:
                    continue
                name = fields[-1]
                if fields[0].startswith("d"):
                    directories.append(name)
                else:
                    files.append(name)

        return sorted(directories), sorted(files)

    def listfiles(self, path: str = "cwd") -> list[str]:
        """List files in an FTP directory.

        :param path: Absolute directory path or ``cwd``.
        :return: Sorted file names.
        """
        if path == "cwd":
            path = f"{self.ftp.pwd().rstrip('/')}/"
        _, files = self.walk(path)
        return files

    def listdirectories(self, path: str = "cwd") -> list[str]:
        """List subdirectories in an FTP directory.

        :param path: Absolute directory path or ``cwd``.
        :return: Sorted directory names.
        """
        if path == "cwd":
            path = f"{self.ftp.pwd().rstrip('/')}/"
        directories, _ = self.walk(path)
        return directories

    def _https_url(self, remote_path: str) -> str:
        """Convert an NCBI FTP URL or path into its HTTPS equivalent."""
        parsed_path = urlparse(remote_path).path
        normalized_path = f"/{parsed_path.lstrip('/')}"
        return f"https://{self._ncbi_host}{normalized_path}"

    def _request(self, remote_path: str, method: str = "GET") -> Request:
        """Build an identified HTTPS request for an NCBI resource."""
        return Request(
            self._https_url(remote_path),
            headers={
                "User-Agent": f"OrthoEvol NCBI downloader ({self.email})",
            },
            method=method,
        )

    def _read_remote_text(self, remote_path: str) -> str:
        """Read a small UTF-8 NCBI resource over HTTPS."""
        with urlopen(
            self._request(remote_path),
            timeout=self._timeout,
        ) as response:
            return response.read().decode("utf-8")

    def _remote_headers(self, remote_path: str) -> dict[str, str]:
        """Return normalized HTTP headers for a remote file."""
        with urlopen(
            self._request(remote_path, method="HEAD"),
            timeout=self._timeout,
        ) as response:
            return {key.lower(): value for key, value in response.headers.items()}

    @staticmethod
    def _md5(file_path: Path) -> str:
        """Calculate an MD5 digest for validation against NCBI sidecars."""
        digest = hashlib.md5()
        with file_path.open("rb") as input_file:
            while chunk := input_file.read(NcbiFTPClient._transfer_chunk_size):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _parse_md5(checksum_text: str) -> str:
        """Extract an MD5 digest from an NCBI checksum file."""
        match = re.search(r"\b[0-9a-fA-F]{32}\b", checksum_text)
        if match is None:
            raise ValueError("NCBI checksum file does not contain an MD5 digest.")
        return match.group(0).lower()

    @staticmethod
    def _marker_matches(marker_path: Path, expected_md5: str) -> bool:
        """Check whether an installed archive marker matches NCBI."""
        if not marker_path.is_file():
            return False
        try:
            observed_md5 = NcbiFTPClient._parse_md5(
                marker_path.read_text(encoding="utf-8")
            )
        except (OSError, ValueError):
            return False
        return observed_md5 == expected_md5

    @staticmethod
    def _write_text_atomic(destination: Path, text: str) -> None:
        """Replace a small text file only after its new content is complete."""
        temporary_path = destination.with_name(f"{destination.name}.part")
        temporary_path.write_text(text, encoding="utf-8")
        temporary_path.replace(destination)

    @staticmethod
    def _local_file_is_current(
        destination: Path,
        remote_headers: dict[str, str],
    ) -> bool:
        """Compare local size and modification time with an HTTPS resource."""
        if not destination.is_file():
            return False

        content_length = remote_headers.get("content-length")
        if content_length is not None and destination.stat().st_size != int(
            content_length
        ):
            return False

        last_modified = remote_headers.get("last-modified")
        if last_modified is None:
            return content_length is not None

        remote_timestamp = parsedate_to_datetime(last_modified).timestamp()
        return destination.stat().st_mtime >= remote_timestamp

    def _download_url(
        self,
        remote_path: str,
        destination: Path,
        expected_md5: str | None = None,
    ) -> Path:
        """Download one file atomically, skipping a current local copy."""
        destination.parent.mkdir(parents=True, exist_ok=True)
        if (
            expected_md5 is not None
            and destination.is_file()
            and self._md5(destination) == expected_md5
        ):
            self.ncbiftp_log.info("%s is current.", destination.name)
            return destination

        remote_headers = self._remote_headers(remote_path)
        if expected_md5 is None and self._local_file_is_current(
            destination,
            remote_headers,
        ):
            self.ncbiftp_log.info("%s is current.", destination.name)
            return destination

        temporary_path = destination.with_name(f"{destination.name}.part")
        digest = hashlib.md5() if expected_md5 is not None else None
        bytes_written = 0

        try:
            with urlopen(
                self._request(remote_path),
                timeout=self._timeout,
            ) as response:
                with temporary_path.open("wb") as output_file:
                    while chunk := response.read(self._transfer_chunk_size):
                        output_file.write(chunk)
                        bytes_written += len(chunk)
                        if digest is not None:
                            digest.update(chunk)

            content_length = remote_headers.get("content-length")
            if content_length is not None and bytes_written != int(content_length):
                raise OSError(
                    f"Incomplete download for {destination.name}: "
                    f"expected {content_length} bytes, received {bytes_written}."
                )
            if digest is not None and digest.hexdigest() != expected_md5:
                raise OSError(f"MD5 validation failed for {destination.name}.")

            temporary_path.replace(destination)
            last_modified = remote_headers.get("last-modified")
            if last_modified is not None:
                remote_timestamp = parsedate_to_datetime(last_modified).timestamp()
                os.utime(destination, (remote_timestamp, remote_timestamp))
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

        self.ncbiftp_log.info("%s was downloaded.", destination.name)
        return destination

    def download_file(
        self,
        filename: str,
        download_path: str | Path = ".",
    ) -> Path:
        """Download a file from the current or specified NCBI directory.

        :param filename: Remote file name or absolute remote path.
        :param download_path: Local destination directory.
        :return: Local file path.
        """
        if filename.startswith("/"):
            remote_path = filename
        else:
            remote_path = f"{self.ftp.pwd().rstrip('/')}/{filename}"
        destination = Path(download_path) / PurePosixPath(filename).name
        return self._download_url(remote_path, destination)

    def _download_pool(
        self,
        remote_paths: Sequence[str],
        download_path: str | Path,
        expected_checksums: dict[str, str] | None = None,
    ) -> list[Path]:
        """Download independent files concurrently over HTTPS."""
        if not remote_paths:
            return []

        destination = Path(download_path)
        expected_checksums = expected_checksums or {}

        def download_remote_file(remote_path: str) -> Path:
            local_path = destination / PurePosixPath(remote_path).name
            return self._download_url(
                remote_path,
                local_path,
                expected_checksums.get(remote_path),
            )

        start_time = perf_counter()
        worker_count = min(self.cpus, len(remote_paths))
        with ThreadPoolExecutor(max_workers=worker_count) as download_pool:
            downloaded_files = list(
                tqdm(
                    download_pool.map(download_remote_file, remote_paths),
                    total=len(remote_paths),
                    unit="file",
                )
            )

        minutes = round((perf_counter() - start_time) / 60, 2)
        self.ncbiftp_log.info(
            "Downloaded %s files in %s minutes.",
            len(downloaded_files),
            minutes,
        )
        return downloaded_files

    @staticmethod
    def _validate_tar_members(
        archive: tarfile.TarFile,
        destination: Path,
    ) -> None:
        """Reject archive members that could write outside the destination."""
        destination_root = destination.resolve()
        for member in archive.getmembers():
            member_path = (destination / member.name).resolve()
            if not member_path.is_relative_to(destination_root):
                raise ValueError(
                    f"Unsafe path {member.name!r} in {archive.name!r}."
                )

            if member.issym() or member.islnk():
                link_path = (member_path.parent / member.linkname).resolve()
                if not link_path.is_relative_to(destination_root):
                    raise ValueError(
                        f"Unsafe link {member.name!r} in {archive.name!r}."
                    )

    def extract_file(
        self,
        file2extract: str | Path,
        download_path: str | Path | None = None,
        remove_archive: bool = True,
    ) -> Path:
        """Extract a gzip or tar-gzip archive into a controlled directory."""
        archive_path = Path(file2extract)
        destination = (
            Path(download_path) if download_path is not None else archive_path.parent
        )
        destination.mkdir(parents=True, exist_ok=True)

        if archive_path.name.endswith(".tar.gz"):
            with tarfile.open(archive_path, mode="r:gz") as archive:
                self._validate_tar_members(archive, destination)
                archive.extractall(destination)
            output_path = destination
        elif archive_path.suffix == ".gz":
            output_path = destination / archive_path.stem
            temporary_path = output_path.with_name(f"{output_path.name}.part")
            try:
                with gzip.open(archive_path, "rb") as compressed_file:
                    with temporary_path.open("wb") as output_file:
                        shutil.copyfileobj(compressed_file, output_file)
                temporary_path.replace(output_path)
            except Exception:
                temporary_path.unlink(missing_ok=True)
                raise
        else:
            raise ValueError(f"Unsupported archive format: {archive_path}")

        if remove_archive:
            archive_path.unlink()
        self.ncbiftp_log.info("%s was extracted.", archive_path.name)
        return output_path

    def _extract_pool(
        self,
        files: Sequence[Path],
        download_path: str | Path | None = None,
    ) -> list[Path]:
        """Extract independent archives concurrently."""
        if not files:
            return []

        start_time = perf_counter()
        worker_count = min(self.cpus, len(files))

        def extract_archive(file_path: Path) -> Path:
            return self.extract_file(file_path, download_path=download_path)

        with ThreadPoolExecutor(max_workers=worker_count) as extract_pool:
            extracted_files = list(
                tqdm(
                    extract_pool.map(extract_archive, files),
                    total=len(files),
                    unit="file",
                )
            )

        minutes = round((perf_counter() - start_time) / 60, 2)
        self.ncbiftp_log.info(
            "Extracted %s files in %s minutes.",
            len(extracted_files),
            minutes,
        )
        return extracted_files

    def _load_blast_metadata(self) -> list[dict[str, Any]]:
        """Load NCBI's authoritative BLAST database manifest."""
        metadata = json.loads(self._read_remote_text(self._blast_metadata_path))
        if not isinstance(metadata, list) or not all(
            isinstance(entry, dict) for entry in metadata
        ):
            raise ValueError("NCBI BLAST metadata has an unexpected structure.")
        return metadata

    @staticmethod
    def _metadata_entry(
        metadata: Sequence[dict[str, Any]],
        database_name: str,
    ) -> dict[str, Any]:
        """Select one exact database entry from the BLAST manifest."""
        for entry in metadata:
            if entry.get("dbname") == database_name:
                return entry
        raise FileNotFoundError(
            f"{database_name!r} is not present in NCBI's BLAST metadata."
        )

    def _blast_archive_paths(
        self,
        database_name: str,
        include_taxonomy: bool = True,
    ) -> list[str]:
        """Resolve every archive required for a version 5 BLAST database."""
        metadata = self._load_blast_metadata()
        database_names = [database_name]
        if include_taxonomy and database_name != "taxdb":
            database_names.append("taxdb")

        remote_paths: list[str] = []
        for selected_name in database_names:
            entry = self._metadata_entry(metadata, selected_name)
            files = entry.get("files")
            if not isinstance(files, list) or not all(
                isinstance(file_url, str) for file_url in files
            ):
                raise ValueError(
                    f"NCBI metadata for {selected_name!r} has no valid files list."
                )
            remote_paths.extend(urlparse(file_url).path for file_url in files)
        return remote_paths

    def _legacy_blast_archive_paths(
        self,
        database_name: str,
        include_taxonomy: bool = True,
    ) -> list[str]:
        """Resolve exact version 4 archive names from its legacy directory."""
        files = self.listfiles(self.blastdbv4_path)
        pattern = re.compile(
            rf"^{re.escape(database_name)}(?:\.\d+)?\.tar\.gz$"
        )
        matching_files = [
            file_name for file_name in files if pattern.fullmatch(file_name)
        ]
        if not matching_files:
            raise FileNotFoundError(
                f"{database_name!r} is not present in NCBI's BLAST v4 directory."
            )

        remote_paths = [
            f"{self.blastdbv4_path}{file_name}" for file_name in matching_files
        ]
        if include_taxonomy and database_name != "taxdb":
            remote_paths.append(f"{self.blastdb_path}taxdb.tar.gz")
        return remote_paths

    def _remote_checksums(
        self,
        remote_paths: Sequence[str],
    ) -> tuple[dict[str, str], dict[str, str]]:
        """Fetch checksum sidecars for a set of NCBI archives."""
        def read_checksum(remote_path: str) -> tuple[str, str, str]:
            checksum_text = self._read_remote_text(f"{remote_path}.md5")
            return remote_path, self._parse_md5(checksum_text), checksum_text

        worker_count = min(self.cpus, len(remote_paths))
        with ThreadPoolExecutor(max_workers=worker_count) as checksum_pool:
            checksum_results = list(checksum_pool.map(read_checksum, remote_paths))

        expected_checksums = {
            remote_path: checksum
            for remote_path, checksum, _ in checksum_results
        }
        checksum_texts = {
            remote_path: checksum_text
            for remote_path, _, checksum_text in checksum_results
        }
        return expected_checksums, checksum_texts

    def getwindowmaskerfiles(
        self,
        taxonomy_ids: Sequence[int | str],
        download_path: str | Path,
    ) -> NoReturn:
        """Reject downloads from NCBI's retired WindowMasker file service."""
        raise OrthoEvolDeprecationWarning(
            "WindowMasker downloads are no longer supported."
        )

    def getblastdb(
        self,
        database_name: str,
        download_path: str | Path,
        v5: bool = True,
        extract: bool = True,
        include_taxonomy: bool = True,
    ) -> None:
        """Download a complete preformatted NCBI BLAST database.

        Version 5 selections come from NCBI's metadata manifest, which avoids
        partial or substring-based database matches.
        """
        if database_name.startswith("est"):
            raise NotImplementedError("EST databases are not supported.")

        destination = Path(download_path)
        destination.mkdir(parents=True, exist_ok=True)
        if v5:
            remote_paths = self._blast_archive_paths(
                database_name,
                include_taxonomy=include_taxonomy,
            )
        else:
            remote_paths = self._legacy_blast_archive_paths(
                database_name,
                include_taxonomy=include_taxonomy,
            )

        expected_checksums, checksum_texts = self._remote_checksums(remote_paths)
        self.files2download = [
            file_name
            for remote_path in remote_paths
            for file_name in (
                PurePosixPath(remote_path).name,
                f"{PurePosixPath(remote_path).name}.md5",
            )
        ]

        paths_to_download: list[str] = []
        for remote_path in remote_paths:
            archive_name = PurePosixPath(remote_path).name
            archive_path = destination / archive_name
            marker_path = destination / f"{archive_name}.md5"
            if (
                extract
                and not archive_path.exists()
                and self._marker_matches(
                    marker_path,
                    expected_checksums[remote_path],
                )
            ):
                self.ncbiftp_log.info("%s is already installed.", archive_name)
                continue
            paths_to_download.append(remote_path)

        downloaded_paths = self._download_pool(
            paths_to_download,
            destination,
            expected_checksums=expected_checksums,
        )

        if extract:
            for remote_path, archive_path in zip(
                paths_to_download,
                downloaded_paths,
                strict=True,
            ):
                self.extract_file(archive_path, download_path=destination)
                marker_path = destination / f"{archive_path.name}.md5"
                self._write_text_atomic(marker_path, checksum_texts[remote_path])
        else:
            for remote_path, archive_path in zip(
                paths_to_download,
                downloaded_paths,
                strict=True,
            ):
                marker_path = destination / f"{archive_path.name}.md5"
                self._write_text_atomic(marker_path, checksum_texts[remote_path])

    def getblastfasta(
        self,
        database_name: str,
        download_path: str | Path,
        extract: bool = True,
    ) -> None:
        """Download one of NCBI's convenience BLAST FASTA archives."""
        if database_name.startswith("est"):
            raise NotImplementedError("EST databases are not supported.")

        archive_name = f"{database_name}.gz"
        if archive_name not in self.listfiles(self.blastfasta_path):
            raise FileNotFoundError(
                f"{database_name!r} is not present in NCBI's BLAST FASTA directory."
            )

        destination = Path(download_path)
        destination.mkdir(parents=True, exist_ok=True)
        remote_path = f"{self.blastfasta_path}{archive_name}"
        expected_checksums, checksum_texts = self._remote_checksums([remote_path])
        archive_path = destination / archive_name
        marker_path = destination / f"{archive_name}.md5"
        extracted_path = destination / database_name
        self.files2download = [archive_name, f"{archive_name}.md5"]

        if (
            extract
            and not archive_path.exists()
            and extracted_path.exists()
            and self._marker_matches(marker_path, expected_checksums[remote_path])
        ):
            self.ncbiftp_log.info("%s is already installed.", archive_name)
            return

        self._download_url(
            remote_path,
            archive_path,
            expected_md5=expected_checksums[remote_path],
        )
        if extract:
            self.extract_file(archive_path, download_path=destination)
        self._write_text_atomic(marker_path, checksum_texts[remote_path])

    def _refseq_marker_path(
        self,
        destination: Path,
        collection_subset: str,
        seqtype: str,
        seqformat: str,
    ) -> Path:
        """Build a stable local marker name for one RefSeq selection."""
        marker_name = (
            f".orthoevol-refseq-{collection_subset}-{seqtype}-{seqformat}.json"
        )
        return destination / marker_name

    def _refseq_checksums(
        self,
        release_number: str,
        file_names: Sequence[str],
    ) -> dict[str, str]:
        """Read checksums for selected files from the RefSeq release catalog."""
        catalog_path = (
            f"{self.refseqrelease_path}release-catalog/"
            f"release{release_number}.files.installed"
        )
        selected_files = set(file_names)
        checksums = {
            file_name: self._parse_md5(checksum)
            for checksum, file_name in (
                line.split("\t", maxsplit=1)
                for line in self._read_remote_text(catalog_path).splitlines()
                if "\t" in line
            )
            if file_name in selected_files
        }
        missing_files = selected_files.difference(checksums)
        if missing_files:
            missing = ", ".join(sorted(missing_files))
            raise ValueError(f"RefSeq release catalog lacks checksums for: {missing}.")
        return checksums

    @staticmethod
    def _refseq_selection_is_current(
        marker_path: Path,
        release_number: str,
        file_names: Sequence[str],
        destination: Path,
        extract: bool,
        checksums: dict[str, str],
    ) -> bool:
        """Check a completed RefSeq selection against the current release."""
        if not marker_path.is_file():
            return False
        try:
            marker_data = json.loads(marker_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return False

        expected_files = [
            str(destination / Path(file_name).with_suffix(""))
            if extract
            else str(destination / file_name)
            for file_name in file_names
        ]
        return (
            marker_data.get("release_number") == release_number
            and marker_data.get("files") == list(file_names)
            and marker_data.get("extract") is extract
            and marker_data.get("checksums") == checksums
            and all(Path(file_path).exists() for file_path in expected_files)
            and (
                extract
                or all(
                    NcbiFTPClient._md5(destination / file_name)
                    == checksums[file_name]
                    for file_name in file_names
                )
            )
        )

    def getrefseqrelease(
        self,
        collection_subset: str,
        seqtype: str,
        seqformat: str,
        download_path: str | Path,
        extract: bool = True,
    ) -> None:
        """Download one molecule and format selection from a RefSeq release."""
        collection_directories = self.listdirectories(self.refseqrelease_path)
        if collection_subset not in collection_directories:
            raise FileNotFoundError(
                f"{collection_subset!r} is not an NCBI RefSeq release collection."
            )

        remote_directory = f"{self.refseqrelease_path}{collection_subset}/"
        release_files = self.listfiles(remote_directory)
        increment_pattern = r"(?:\d+(?:\.\d+)?|wp_protein\.\d+)"
        pattern = re.compile(
            rf"^{re.escape(collection_subset)}\.{increment_pattern}\."
            rf"{re.escape(seqtype)}\.{re.escape(seqformat)}\.gz$"
        )
        selected_files = [
            file_name
            for file_name in release_files
            if pattern.fullmatch(file_name)
        ]
        if not selected_files:
            raise FileNotFoundError(
                "No RefSeq release files matched "
                f"{collection_subset=}, {seqtype=}, {seqformat=}."
            )

        destination = Path(download_path)
        destination.mkdir(parents=True, exist_ok=True)
        self.files2download = selected_files
        release_number = self._read_remote_text(
            self.refseq_release_number_path
        ).strip()
        checksums = self._refseq_checksums(release_number, selected_files)
        marker_path = self._refseq_marker_path(
            destination,
            collection_subset,
            seqtype,
            seqformat,
        )
        if self._refseq_selection_is_current(
            marker_path,
            release_number,
            selected_files,
            destination,
            extract,
            checksums,
        ):
            self.ncbiftp_log.info(
                "RefSeq release %s selection is current.",
                release_number,
            )
            return

        remote_paths = [
            f"{remote_directory}{file_name}" for file_name in selected_files
        ]
        expected_checksums = {
            remote_path: checksums[PurePosixPath(remote_path).name]
            for remote_path in remote_paths
        }
        downloaded_files = self._download_pool(
            remote_paths,
            destination,
            expected_checksums,
        )
        if extract:
            self._extract_pool(downloaded_files, download_path=destination)

        marker_data = {
            "release_number": release_number,
            "files": selected_files,
            "extract": extract,
            "checksums": checksums,
        }
        self._write_text_atomic(
            marker_path,
            f"{json.dumps(marker_data, indent=2, sort_keys=True)}\n",
        )

    def updatedb(
        self,
        database_path: str | Path | None = None,
        update_days: int = 7,
    ) -> None:
        """Refresh recognized local databases after a minimum age."""
        if update_days < 0:
            raise ValueError("update_days cannot be negative.")

        destination = Path.cwd() if database_path is None else Path(database_path)
        blast_aliases = sorted(
            {
                file_path.stem
                for suffix in ("*.nal", "*.pal")
                for file_path in destination.glob(suffix)
            }
        )
        refseq_files = sorted(destination.glob("*.gbff"))
        recognized_files = [
            *(
                file_path
                for suffix in ("*.nal", "*.pal")
                for file_path in destination.glob(suffix)
            ),
            *refseq_files,
        ]
        if not recognized_files:
            raise NotImplementedError(
                "No BLAST alias or extracted RefSeq files were found."
            )

        newest_mtime = max(file_path.stat().st_mtime for file_path in recognized_files)
        time_elapsed = datetime.now() - datetime.fromtimestamp(newest_mtime)
        if time_elapsed.days < update_days:
            self.ncbiftp_log.info(
                "The local database is still within its update window."
            )
            return

        if blast_aliases:
            for database_name in blast_aliases:
                self.getblastdb(database_name, destination, extract=True)
            return

        file_parts = refseq_files[0].name.split(".")
        if len(file_parts) < 4:
            raise ValueError(
                f"Cannot infer RefSeq selection from {refseq_files[0].name!r}."
            )
        self.getrefseqrelease(
            collection_subset=file_parts[0],
            seqtype=file_parts[-2],
            seqformat=file_parts[-1],
            download_path=destination,
            extract=True,
        )
