"""Explicit preparation and offline access to ETE taxonomy databases."""

import logging
import sqlite3
import sys
from collections.abc import Iterator, Sequence
from contextlib import closing, contextmanager
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from urllib.request import urlretrieve

from ete3 import NCBITaxa
from ete3.ncbi_taxonomy.ncbiquery import DB_VERSION

logger = logging.getLogger(__name__)
TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


def validate_taxonomy_database(dbfile: str | Path | None = None) -> Path:
    """Require a populated, compatible database without creating or updating it."""
    path = Path(dbfile or "~/.etetoolkit/taxa.sqlite").expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Taxonomy database not found: {path}. "
            "Run prepare_taxonomy_database() before starting the analysis."
        )
    try:
        with closing(sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True)) as db:
            version = db.execute("SELECT version FROM stats").fetchone()
            root = db.execute("SELECT track FROM species WHERE taxid = 1").fetchone()
            populated = db.execute(
                "SELECT taxid FROM species WHERE taxid != 1 LIMIT 1"
            ).fetchone()
            if version != (DB_VERSION,) or not root or root[0] != "1" or not populated:
                raise ValueError("Missing taxonomy rows or incompatible ETE schema")
    except (sqlite3.DatabaseError, ValueError) as error:
        raise ValueError(
            f"Invalid taxonomy database: {path}. "
            "Rebuild it with prepare_taxonomy_database(..., overwrite=True)."
        ) from error
    return path


@contextmanager
def taxonomy_database(dbfile: str | Path | None = None) -> Iterator[NCBITaxa]:
    """Open a validated database for offline lookups and close its connection."""
    path = validate_taxonomy_database(dbfile)
    taxonomy = NCBITaxa(dbfile=str(path), update=False)
    try:
        yield taxonomy
    finally:
        taxonomy.db.close()


def resolve_species(taxonomy: NCBITaxa, species: Sequence[str]) -> dict[str, int]:
    """Resolve every input species; missing or ambiguous names are errors."""
    names = [name.replace("_", " ") for name in species]
    matches = taxonomy.get_name_translator(names)
    unresolved = [name for name in names if len(matches.get(name, [])) != 1]
    if unresolved:
        raise ValueError(
            f"Missing or ambiguous taxonomy names: {', '.join(unresolved)}"
        )
    return {original: matches[name][0] for original, name in zip(species, names)}


def prepare_taxonomy_database(
    dbfile: str | Path,
    *,
    taxdump_file: str | Path | None = None,
    required_species: Sequence[str] = (),
    overwrite: bool = False,
) -> Path:
    """Build and validate taxonomy before publishing it at the requested path.

    Reuse an existing valid database unless overwrite=True. A local NCBI taxdump
    can be supplied for reproducibility; otherwise download it over HTTPS.
    Run preparation before analyses, without concurrent writers/readers: ETE's
    SQLite database and traversal cache cannot be replaced as one transaction.
    Failed builds leave the existing database and cache untouched.
    """
    destination = Path(dbfile).expanduser().resolve()
    if destination.exists() and not overwrite:
        with taxonomy_database(destination) as taxonomy:
            resolve_species(taxonomy, required_species)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".taxonomy-", dir=destination.parent) as directory:
        work = Path(directory)
        archive = (
            Path(taxdump_file).expanduser().resolve()
            if taxdump_file
            else work / "taxdump.tar.gz"
        )
        if taxdump_file is None:
            logger.info("Downloading NCBI taxonomy to %s", archive)
            urlretrieve(TAXDUMP_URL, archive)
        if not archive.is_file():
            raise FileNotFoundError(archive)
        candidate = work / "taxa.sqlite"
        # ETE writes fixed temporary filenames in cwd. Isolate those files without
        # changing the calling process's working directory or shared state.
        run(
            [
                sys.executable,
                "-c",
                "import sys; from ete3 import NCBITaxa; "
                "ncbi = NCBITaxa(dbfile=sys.argv[1], taxdump_file=sys.argv[2]); "
                "ncbi.db.close()",
                str(candidate),
                str(archive),
            ],
            cwd=work,
            check=True,
        )
        with taxonomy_database(candidate) as taxonomy:
            if taxonomy.db.execute("PRAGMA quick_check").fetchone() != ("ok",):
                raise ValueError("Taxonomy database failed SQLite integrity validation")
            resolve_species(taxonomy, required_species)
        Path(f"{candidate}.traverse.pkl").replace(f"{destination}.traverse.pkl")
        candidate.replace(destination)
    return destination
