"""Validated models for OrthoEvolution YAML configuration files."""

from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar, Literal

import yaml
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StrictBool,
    StrictStr,
)


class ConfigurationFileError(ValueError):
    """Report configuration files that cannot be parsed structurally."""


class ConfigurationSection(BaseModel):
    """Base model that rejects unknown public configuration fields."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class ManagementConfig(ConfigurationSection):
    """Validate inputs used to construct ``ProjectManagement``."""

    repo: StrictStr | None
    user: StrictStr | None
    project: StrictStr
    research: StrictStr | None = None
    research_type: StrictStr | None = None
    app: StrictStr | None = None
    home: Path | None = None
    new_repo: StrictBool | None = None
    new_user: StrictBool | None = None
    new_project: StrictBool | None = None
    new_research: StrictBool | None = None
    new_app: StrictBool | None = None
    new_website: StrictBool | None = None
    new_db: StrictBool | None = None


class StrategyFlagsConfig(ConfigurationSection):
    """Validate execution flags inherited by nested database strategies."""

    configure_flag: StrictBool | None = None
    archive_flag: StrictBool | None = None
    delete_flag: StrictBool | None = None


class DatabaseStrategyConfig(StrategyFlagsConfig):
    """Validate flags and paths shared by concrete database strategies."""

    database_path: Path | None = None
    archive_path: Path | None = None


class RefSeqReleaseStrategyConfig(DatabaseStrategyConfig):
    """Validate active inputs for an NCBI RefSeq release selection."""

    collection_subset: StrictStr | None = None
    seqtype: StrictStr | None = None
    seqformat: StrictStr | None = None
    template_flag: StrictBool | None = None
    download_flag: StrictBool | None = None


class NcbiBlastStrategyConfig(DatabaseStrategyConfig):
    """Validate the BLAST strategy and its required database child."""

    ncbi_blast_db: DatabaseStrategyConfig = Field(alias="NCBI_blast_db")


class NcbiStrategyConfig(DatabaseStrategyConfig):
    """Validate the complete set of children required by the NCBI strategy."""

    ncbi_blast: NcbiBlastStrategyConfig = Field(alias="NCBI_blast")
    ncbi_pub_taxonomy: DatabaseStrategyConfig = Field(
        alias="NCBI_pub_taxonomy",
    )
    ncbi_refseq_release: RefSeqReleaseStrategyConfig = Field(
        alias="NCBI_refseq_release",
    )


class ItisStrategyConfig(DatabaseStrategyConfig):
    """Validate the taxonomy child required by the ITIS strategy."""

    itis_taxonomy: DatabaseStrategyConfig = Field(alias="ITIS_taxonomy")


class FullDatabaseStrategyConfig(StrategyFlagsConfig):
    """Validate both children required by the full database strategy."""

    ncbi: NcbiStrategyConfig = Field(alias="NCBI")
    itis: ItisStrategyConfig = Field(alias="ITIS")


class DatabaseConfig(ConfigurationSection):
    """Validate base database inputs and recognized strategy names."""

    email: StrictStr
    driver: StrictStr
    project: StrictStr | None = None
    project_path: Path | None = None
    blast: StrictBool | None = None
    ftp_flag: StrictBool | None = None
    full: FullDatabaseStrategyConfig | None = Field(default=None, alias="Full")
    ncbi: NcbiStrategyConfig | None = Field(default=None, alias="NCBI")
    ncbi_blast: NcbiBlastStrategyConfig | None = Field(
        default=None,
        alias="NCBI_blast",
    )
    ncbi_blast_db: DatabaseStrategyConfig | None = Field(
        default=None,
        alias="NCBI_blast_db",
    )
    ncbi_pub_taxonomy: DatabaseStrategyConfig | None = Field(
        default=None,
        alias="NCBI_pub_taxonomy",
    )
    ncbi_refseq_release: RefSeqReleaseStrategyConfig | None = Field(
        default=None,
        alias="NCBI_refseq_release",
    )
    itis: ItisStrategyConfig | None = Field(default=None, alias="ITIS")
    itis_taxonomy: DatabaseStrategyConfig | None = Field(
        default=None,
        alias="ITIS_taxonomy",
    )


class ComparativeGeneticsConfig(ConfigurationSection):
    """Validate established comparative-genetics configuration values."""

    project: StrictStr | None = None
    project_path: Path | None = None
    acc_file: StrictStr | None = None
    taxon_file: Path | None = None
    ref_species: StrictStr | None = None
    pre_blast: StrictBool | None = None
    post_blast: StrictBool | None = None
    hgnc: StrictBool | StrictStr | Path | None = None
    copy_from_package: StrictBool | None = None
    save_data: StrictBool | None = None
    template: Any = None
    go_list: Any = None
    maf: StrictStr | None = Field(default=None, alias="MAF")


class BlastConfig(ConfigurationSection):
    """Validate established nucleotide-BLAST configuration values."""

    project: StrictStr | None = None
    project_path: Path | None = None
    method: Literal[1, 2] | None = None
    template: Any = None
    save_data: StrictBool | None = None
    acc_file: StrictStr | None = None
    copy_from_package: StrictBool | None = None
    auto_start: StrictBool | None = None


class GenBankConfig(ConfigurationSection):
    """Validate established GenBank configuration flags."""

    project: StrictStr | None = None
    project_path: Path | None = None
    solo: StrictBool | None = None
    multi: StrictBool | None = None
    archive: StrictBool | None = None
    min_fasta: StrictBool | None = None


class AlignmentConfig(ConfigurationSection):
    """Validate alignment selection while retaining tool-specific options."""

    aln_program: StrictStr | None = None
    guidance_config: dict[str, Any] | None = Field(
        default=None,
        alias="Guidance_config",
    )
    pal2nal_config: dict[str, Any] | None = Field(
        default=None,
        alias="Pal2Nal_config",
    )
    clustalo_config: dict[str, Any] | None = Field(
        default=None,
        alias="ClustalO_config",
    )
    legacy_clustalo_config: dict[str, Any] | None = Field(
        default=None,
        alias="Clustalo_config",
    )


class PipelineConfig(BaseModel):
    """Represent the supported top-level OrthoEvolution configuration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    section_names: ClassVar[tuple[str, ...]] = (
        "Management_config",
        "Database_config",
        "CompGenAnalysis_config",
        "BLASTn_config",
        "GenBank_config",
        "Alignment_config",
    )

    management: ManagementConfig | None = Field(
        default=None,
        alias="Management_config",
    )
    database: DatabaseConfig | None = Field(
        default=None,
        alias="Database_config",
    )
    comparative_genetics: ComparativeGeneticsConfig | None = Field(
        default=None,
        alias="CompGenAnalysis_config",
    )
    blastn: BlastConfig | None = Field(default=None, alias="BLASTn_config")
    genbank: GenBankConfig | None = Field(default=None, alias="GenBank_config")
    alignment: AlignmentConfig | None = Field(
        default=None,
        alias="Alignment_config",
    )

    def as_legacy_dict(self) -> dict[str, dict[str, Any]]:
        """Return validated values using the keys expected by legacy classes."""
        return self.model_dump(
            by_alias=True,
            exclude_none=True,
            exclude_unset=True,
        )


def load_pipeline_config(config_file: str | Path) -> PipelineConfig:
    """Load and validate an OrthoEvolution YAML configuration file."""
    config_path = Path(config_file)
    try:
        raw_config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise ConfigurationFileError(
            f"Unable to parse configuration file {config_path}: {error}"
        ) from error

    if raw_config is None:
        raise ConfigurationFileError(f"Configuration file {config_path} is empty.")
    if not isinstance(raw_config, Mapping):
        raise ConfigurationFileError(
            f"Configuration file {config_path} must contain a YAML mapping."
        )

    return PipelineConfig.model_validate(raw_config)
