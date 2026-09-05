"""Focused tests for the SQLite BioSQL workflow."""

from pathlib import Path
from unittest import mock

import pytest

from OrthoEvol.Manager.biosql.biosql import BaseBioSQL, SQLiteBioSQL


def build_sqlite(tmp_path: Path) -> SQLiteBioSQL:
    """Build only the paths and collaborators used by SQLite methods."""
    biosql = object.__new__(SQLiteBioSQL)
    biosql.template_abs_path = tmp_path / "template.db"
    biosql.databases_path = tmp_path / "databases"
    biosql.database_name = Path("project.db")
    biosql.driver = "sqlite3"
    biosql.schema_cmd = "sqlite3 %s -echo"
    biosql.schema_file = "biosqldb-sqlite.sql"
    biosql.taxon_cmd = "%s --dbname %s --driver %s --download false --directory %s"
    biosql.ncbi_taxon_script = tmp_path / "load_taxonomy.pl"
    biosql.biosqllog = mock.Mock()
    biosql.configure_new_database = mock.Mock()
    return biosql


def test_configure_new_database_adds_schema_redirection() -> None:
    biosql = object.__new__(BaseBioSQL)
    biosql.biosql_proc = mock.Mock()

    biosql.configure_new_database("sqlite3 template.db", "schema.sql")

    biosql.biosql_proc.assert_called_once_with(
        cmd="sqlite3 template.db < schema.sql",
        stdout=mock.ANY,
        stderr=mock.ANY,
        shell=True,
    )


def test_load_sqlite_schema_only_for_missing_template(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)

    with mock.patch(
        "OrthoEvol.Manager.biosql.biosql.package_resource_path",
        return_value=tmp_path / "schema.sql",
    ):
        biosql.load_sqlite_schema()
        biosql.template_abs_path.touch()
        biosql.load_sqlite_schema()

    biosql.configure_new_database.assert_called_once_with(
        f"sqlite3 {biosql.template_abs_path} -echo", tmp_path / "schema.sql"
    )


@pytest.mark.parametrize("size", [120_000, 200_000])
def test_load_sqlite_taxonomy_uses_template_size_threshold(
    size: int, tmp_path: Path
) -> None:
    biosql = build_sqlite(tmp_path)
    biosql.template_abs_path.write_bytes(b"0" * size)

    biosql.load_sqlite_taxonomy()

    if size == 120_000:
        biosql.configure_new_database.assert_called_once()
    else:
        biosql.configure_new_database.assert_not_called()


def test_create_template_database_runs_each_stage(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)
    biosql.load_sqlite_schema = mock.Mock()
    biosql.create_executable_scripts = mock.Mock()
    biosql.load_sqlite_taxonomy = mock.Mock()

    biosql.create_template_database()

    biosql.load_sqlite_schema.assert_called_once_with()
    biosql.create_executable_scripts.assert_called_once_with()
    biosql.load_sqlite_taxonomy.assert_called_once_with()


def test_copy_template_database_creates_and_copies_template(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)
    biosql.create_template_database = mock.Mock(
        side_effect=biosql.template_abs_path.touch
    )
    process = mock.Mock()

    with mock.patch(
        "OrthoEvol.Manager.biosql.biosql.sp.Popen", return_value=process
    ) as popen:
        biosql.copy_template_database(tmp_path)

    biosql.create_template_database.assert_called_once_with()
    popen.assert_called_once_with(
        ["cp", str(biosql.template_abs_path), str(tmp_path / "project.db")]
    )
    process.wait.assert_called_once_with()


def test_update_sqlite_taxonomy_uses_existing_template(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)
    biosql.template_abs_path.touch()

    biosql.update_sqlite_taxonomy()

    command = biosql.configure_new_database.call_args.args[0]
    assert str(biosql.template_abs_path) in command
    assert "--driver SQLite" in command


def test_upload_files_loads_and_commits_records(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)
    database = tmp_path / biosql.database_name
    database.touch()
    sub_database = mock.Mock()
    sub_database.load.return_value = 2
    server = mock.MagicMock()
    server.keys.return_value = []
    server.__getitem__.return_value = sub_database

    with (
        mock.patch(
            "OrthoEvol.Manager.biosql.biosql.BioSeqDatabase.open_database",
            return_value=server,
        ),
        mock.patch(
            "OrthoEvol.Manager.biosql.biosql.SeqIO.parse",
            return_value=iter(()),
        ) as parse,
    ):
        biosql.upload_files("rna", "genbank", tmp_path, ["records.gbff"])

    server.new_database.assert_called_once_with("rna")
    parse.assert_called_once_with(str(tmp_path / "records.gbff"), "genbank")
    sub_database.load.assert_called_once()
    server.commit.assert_called_once_with()


def test_create_executable_scripts_only_changes_perl_files(tmp_path: Path) -> None:
    biosql = object.__new__(BaseBioSQL)
    biosql.scripts = tmp_path
    biosql.biosqllog = mock.Mock()
    perl_script = tmp_path / "load_taxonomy.pl"
    text_file = tmp_path / "README.txt"
    perl_script.touch(mode=0o600)
    text_file.touch(mode=0o600)

    biosql.create_executable_scripts()

    assert perl_script.stat().st_mode & 0o777 == 0o755
    assert text_file.stat().st_mode & 0o777 == 0o600


def test_base_biosql_configures_standalone_paths(tmp_path: Path) -> None:
    biosql = BaseBioSQL(
        database_name="project.db",
        template_name="template.db",
        project="project",
        project_path=tmp_path,
        proj_mana=None,
    )

    assert biosql.project_path == tmp_path / "project"
    assert biosql.template_abs_path == tmp_path / "project" / "index" / "template.db"
    assert biosql.databases_path == tmp_path / "databases"


def test_load_sqlite_taxonomy_reports_missing_template(tmp_path: Path) -> None:
    biosql = build_sqlite(tmp_path)

    biosql.load_sqlite_taxonomy()

    biosql.biosqllog.error.assert_called_once()
