"""Explicit preparation and offline access to ETE taxonomy databases."""

import sqlite3
from collections.abc import Iterator
from contextlib import closing, contextmanager
from pathlib import Path

from ete3 import NCBITaxa
from ete3.ncbi_taxonomy.ncbiquery import DB_VERSION

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
