#!/usr/bin/env python
"""Download a current preformatted NCBI BLAST database."""

import argparse
from pathlib import Path

from OrthoEvol.Tools.ftp import NcbiFTPClient


def download_blast_database(
    email: str,
    database_name: str,
    download_path: Path,
    max_workers: int = 8,
) -> Path:
    """Download and extract one database selected from NCBI's manifest."""
    client = NcbiFTPClient(email=email, max_workers=max_workers)
    try:
        return client.getblastdb(
            database_name=database_name,
            download_path=download_path,
        )
    finally:
        # Close the listing connection even when an HTTPS transfer fails.
        client.close_connection()


def parse_arguments() -> argparse.Namespace:
    """Parse inputs for the standalone download example."""
    parser = argparse.ArgumentParser(
        description="Download a current preformatted NCBI BLAST database.",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Contact email used for anonymous NCBI FTP access.",
    )
    parser.add_argument(
        "--database-name",
        required=True,
        help="Exact database name from NCBI's BLAST metadata manifest.",
    )
    parser.add_argument(
        "--download-path",
        type=Path,
        default=Path.cwd(),
        help="Destination directory. Defaults to the current directory.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Maximum number of concurrent downloads.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the standalone BLAST database download."""
    arguments = parse_arguments()
    download_blast_database(
        email=arguments.email,
        database_name=arguments.database_name,
        download_path=arguments.download_path,
        max_workers=arguments.max_workers,
    )


if __name__ == "__main__":
    main()
