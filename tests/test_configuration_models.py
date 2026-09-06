"""Tests for validated OrthoEvolution configuration input."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from OrthoEvol.Manager.config import yml
from OrthoEvol.config import (
    ConfigurationFileError,
    load_pipeline_config,
)
from OrthoEvol.resources import package_resource_path


def write_config(tmp_path: Path, content: str) -> Path:
    """Write a small configuration fixture and return its path."""
    config_file = tmp_path / "config.yml"
    config_file.write_text(content, encoding="utf-8")
    return config_file


def test_loads_valid_configuration(tmp_path: Path) -> None:
    config_file = write_config(
        tmp_path,
        """Management_config:
  repo: repository
  user: researcher
  project: orthology
GenBank_config:
  solo: true
""",
    )

    configuration = load_pipeline_config(config_file).as_legacy_dict()

    assert configuration["Management_config"]["project"] == "orthology"
    assert configuration["GenBank_config"]["solo"] is True


def test_rejects_unknown_section_field(tmp_path: Path) -> None:
    config_file = write_config(
        tmp_path,
        "GenBank_config:\n  unknown_option: true\n",
    )

    with pytest.raises(ValidationError, match="unknown_option"):
        load_pipeline_config(config_file)


def test_rejects_unknown_top_level_section(tmp_path: Path) -> None:
    config_file = write_config(tmp_path, "Unknown_config: {}\n")

    with pytest.raises(ValidationError, match="Unknown_config"):
        load_pipeline_config(config_file)


def test_rejects_non_mapping_section(tmp_path: Path) -> None:
    config_file = write_config(tmp_path, "GenBank_config: enabled\n")

    with pytest.raises(ValidationError, match="GenBank_config"):
        load_pipeline_config(config_file)


@pytest.mark.parametrize(
    ("content", "invalid_field"),
    [
        ("GenBank_config:\n  solo: 'true'\n", "solo"),
        ("BLASTn_config:\n  method: '1'\n", "method"),
        ("BLASTn_config:\n  method: 3\n", "method"),
    ],
)
def test_rejects_coercible_scalar_values(
    tmp_path: Path,
    content: str,
    invalid_field: str,
) -> None:
    config_file = write_config(tmp_path, content)

    with pytest.raises(ValidationError, match=invalid_field):
        load_pipeline_config(config_file)


def test_reports_missing_configuration_file(tmp_path: Path) -> None:
    missing_file = tmp_path / "missing.yml"

    with pytest.raises(FileNotFoundError, match="missing.yml"):
        load_pipeline_config(missing_file)


def test_reports_malformed_yaml(tmp_path: Path) -> None:
    config_file = write_config(tmp_path, "GenBank_config: [\n")

    with pytest.raises(ConfigurationFileError, match="Unable to parse"):
        load_pipeline_config(config_file)


def test_reports_empty_yaml(tmp_path: Path) -> None:
    config_file = write_config(tmp_path, "")

    with pytest.raises(ConfigurationFileError, match="is empty"):
        load_pipeline_config(config_file)


def test_rejects_non_mapping_yaml(tmp_path: Path) -> None:
    config_file = write_config(tmp_path, "- GenBank_config\n")

    with pytest.raises(ConfigurationFileError, match="YAML mapping"):
        load_pipeline_config(config_file)


def test_converts_explicit_database_path(tmp_path: Path) -> None:
    config_file = write_config(
        tmp_path,
        """Database_config:
  email: test@example.com
  driver: sqlite3
  project_path: databases
""",
    )

    database_config = load_pipeline_config(config_file).as_legacy_dict()[
        "Database_config"
    ]

    assert database_config["project_path"] == Path("databases")


def test_rejects_unknown_database_strategy(tmp_path: Path) -> None:
    config_file = write_config(
        tmp_path,
        """Database_config:
  email: test@example.com
  driver: sqlite3
  Unsupported_database: {}
""",
    )

    with pytest.raises(ValidationError, match="Unsupported_database"):
        load_pipeline_config(config_file)


@pytest.mark.parametrize(
    ("nested_configuration", "invalid_field"),
    [
        (
            """Full:
    NCBI:
      Unsupported_database: {}
""",
            "Unsupported_database",
        ),
        (
            """NCBI_blast:
    NCBI_blast_db:
      configure_flag: 'true'
""",
            "configure_flag",
        ),
    ],
)
def test_rejects_invalid_nested_database_strategy_inputs(
    tmp_path: Path,
    nested_configuration: str,
    invalid_field: str,
) -> None:
    config_file = write_config(
        tmp_path,
        "Database_config:\n"
        "  email: test@example.com\n"
        "  driver: sqlite3\n"
        f"  {nested_configuration}",
    )

    with pytest.raises(ValidationError, match=invalid_field):
        load_pipeline_config(config_file)


def test_nested_database_models_preserve_dispatcher_keys(tmp_path: Path) -> None:
    config_file = write_config(
        tmp_path,
        """Database_config:
  email: test@example.com
  driver: sqlite3
  NCBI:
    NCBI_blast:
      NCBI_blast_db:
        configure_flag: false
    NCBI_pub_taxonomy: {}
    NCBI_refseq_release:
      collection_subset: vertebrate_mammalian
      seqtype: rna
      seqformat: gbff
""",
    )

    database_config = load_pipeline_config(config_file).as_legacy_dict()[
        "Database_config"
    ]

    assert database_config["NCBI"]["NCBI_blast"]["NCBI_blast_db"] == {
        "configure_flag": False,
    }
    assert database_config["NCBI"]["NCBI_refseq_release"] == {
        "collection_subset": "vertebrate_mammalian",
        "seqtype": "rna",
        "seqformat": "gbff",
    }


@pytest.mark.parametrize(
    "config_name",
    [
        "pipeline.yml",
        "databases.yml",
        "initialize_new.yml",
        "initialize_old.yml",
    ],
)
def test_packaged_configurations_are_valid(config_name: str) -> None:
    config_file = package_resource_path(yml, config_name)

    assert load_pipeline_config(config_file).as_legacy_dict()
