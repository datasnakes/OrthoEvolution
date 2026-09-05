# NCBI download tools

This package provides a generic FTP connection class and an NCBI-specific
client. `NcbiFTPClient` uses FTP only to inspect directories. File transfers use
HTTPS so concurrent downloads do not share a stateful FTP connection.

The client supports:

- Current version 5 BLAST databases from `/blast/db/`
- Legacy version 4 BLAST databases from `/blast/db/v4/`
- NCBI's convenience BLAST FASTA files
- Bulk RefSeq release selections
- FTP file and directory listings

NCBI no longer publishes WindowMasker files. Calling
`getwindowmaskerfiles()` raises `OrthoEvolDeprecationWarning`.

## BLAST database download

```python
from pathlib import Path

from OrthoEvol.Tools.ftp import NcbiFTPClient


download_path = Path("databases") / "NCBI" / "blast" / "db"
client = NcbiFTPClient(email="researcher@example.org", max_workers=4)

try:
    client.getblastdb(
        database_name="refseq_rna",
        download_path=download_path,
        v5=True,
        extract=True,
    )
finally:
    client.close_connection()
```

Version 5 selections use NCBI's
[`blastdb-metadata-1-1.json`](https://ftp.ncbi.nlm.nih.gov/blast/db/blastdb-metadata-1-1.json)
manifest. Exact database identities are matched, and every volume listed by
NCBI is downloaded. Each archive is verified against its `.md5` sidecar before
extraction. The sidecar remains as the local installation marker.

Set `v5=False` only when a legacy version 4 database is specifically required.

## RefSeq release download

```python
from pathlib import Path

from OrthoEvol.Tools.ftp import NcbiFTPClient


download_path = Path("databases") / "NCBI" / "refseq" / "release"
client = NcbiFTPClient(email="researcher@example.org", max_workers=4)

try:
    client.getrefseqrelease(
        collection_subset="vertebrate_mammalian",
        seqtype="rna",
        seqformat="gbff",
        download_path=download_path,
        extract=True,
    )
finally:
    client.close_connection()
```

RefSeq filenames are matched using NCBI's documented structure:
`collection.increment[.subpart].molecule.format.gz`. The special
`collection.wp_protein.increment.protein.format.gz` form is also included. A
local selection marker records the completed NCBI release number. Repeating the
same request skips a complete selection from the current release. Archives are
verified against NCBI's release-catalog checksums before extraction.

## BLAST FASTA download

```python
client.getblastfasta(
    database_name="swissprot",
    download_path=Path("databases") / "NCBI" / "blast" / "fasta",
)
```

The FASTA directory is a convenience snapshot and can lag the preformatted
BLAST databases. For sequences from the current preformatted database, NCBI
recommends downloading the database and exporting sequences with `blastdbcmd`.

## Directory listings

```python
directories = client.listdirectories("/blast/db/")
files = client.listfiles("/blast/db/")
```

Paths must begin and end with `/`.

## Transfer behavior

- Downloads use temporary `.part` files and replace destinations atomically.
- Existing validated BLAST archives or installation markers are reused.
- Existing RefSeq files are reused when their checksums or extracted outputs agree.
- Archive paths are validated before extraction.
- Network errors and checksum mismatches leave prior complete files intact.

See [NCBIREADME.md](NCBIREADME.md) for the NCBI source-of-truth links and
format notes.
