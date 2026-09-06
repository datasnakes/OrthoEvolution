# Current NCBI download layout

Static database tables become inaccurate as NCBI adds, renames, and retires
resources. This client therefore reads NCBI's live directory listings and BLAST
metadata instead of maintaining a bundled database-name catalog.

## Preformatted BLAST databases

NCBI serves its current BLAST databases from
[`/blast/db/`](https://ftp.ncbi.nlm.nih.gov/blast/db/). OrthoEvolution follows
that current location and does not expose a database-format version switch.

The authoritative manifest is
[`blastdb-metadata-1-1.json`](https://ftp.ncbi.nlm.nih.gov/blast/db/blastdb-metadata-1-1.json).
Each entry supplies:

- The exact database name
- Database type and description
- Every required archive volume
- Update time and compressed size

All volumes are required for a multi-volume database. `getblastdb()` uses the
manifest and also installs `taxdb` by default. Taxonomy data are needed when
BLAST reports translate taxonomy identifiers into names.

NCBI's supported command-line alternative is `update_blastdb.pl`, distributed
with BLAST+. See the
[NCBI BLAST database download guide](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

## BLAST FASTA files

The [`/blast/db/FASTA/`](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
directory contains a small convenience set such as `nr.gz`, `nt.gz`,
`pdbaa.gz`, and `swissprot.gz`. These files are not the preferred source for a
current local BLAST database. NCBI recommends downloading the preformatted
database and using `blastdbcmd` when FASTA output is required.

## RefSeq releases

Bulk RefSeq releases remain available under
[`/refseq/release/`](https://ftp.ncbi.nlm.nih.gov/refseq/release/). The current
release number is published in `RELEASE_NUMBER`.

Most sequence files follow this structure:

```text
collection.increment[.subpart].molecule.format.gz
```

Non-redundant `WP_` protein files follow this special structure:

```text
collection.wp_protein.increment.protein.format.gz
```

Examples include:

```text
vertebrate_mammalian.1.rna.gbff.gz
vertebrate_mammalian.1.1.genomic.fna.gz
vertebrate_mammalian.1.protein.faa.gz
```

Valid molecule values include `genomic`, `rna`, and `protein`. Common formats
include `gbff`, `gpff`, `fna`, and `faa`. The complete definitions and edge
cases are documented in NCBI's
[`RefSeq release README`](https://ftp.ncbi.nlm.nih.gov/refseq/release/README).

For assembly-specific genome downloads, NCBI recommends NCBI Datasets or the
assembly paths in the Genomes FTP hierarchy rather than scanning a complete
RefSeq release.

## Protocols

NCBI continues to support FTP and HTTPS. NCBI retired rsync access on June 1,
2026. This package does not use rsync. See the
[NCBI protocol announcement](https://ncbiinsights.ncbi.nlm.nih.gov/2026/03/25/retire-rsync-support-ftp-downloads/).

## Integrity and restart behavior

BLAST archives are checked against NCBI's MD5 sidecars. Downloads are written to
temporary `.part` files, so a failed transfer does not replace a complete local
file. After successful extraction, the checksum sidecar is retained as the
installation marker. RefSeq selections use NCBI's release number and verify
that every expected output exists before treating the local selection as
current.
