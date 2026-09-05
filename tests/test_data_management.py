"""Unit tests for legacy data-management dispatch."""

from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from OrthoEvol.Manager.data_management import DataMana
from OrthoEvol.Orthologs.Blast.blast import OrthoBlastN


def test_configure_routes_each_enabled_stage(tmp_path: Path) -> None:
    config_file = tmp_path / "pipeline.yml"
    config_file.write_text(
        """Management_config:
  repo: null
  user: null
  project: example
Database_config:
  driver: sqlite
CompGenAnalysis_config:
  taxon_file: taxa.csv
BLASTn_config:
  method: 1
GenBank_config:
  solo: true
Alignment_config:
  Clustalo_config: {}
""",
        encoding="utf-8",
    )
    manager = DataMana()
    manager.database = mock.Mock()
    manager.blast = mock.Mock()
    manager.genbank = mock.Mock()
    manager.align = mock.Mock()
    project_manager = mock.Mock()

    with mock.patch(
        "OrthoEvol.Manager.data_management.ProjectManagement",
        return_value=project_manager,
    ):
        manager.configure(config_file)

    manager.database.assert_called_once_with(
        project_manager,
        {"driver": "sqlite"},
    )
    manager.blast.assert_called_once_with(
        project_manager,
        {"method": 1, "taxon_file": "taxa.csv"},
    )
    manager.genbank.assert_called_once_with(project_manager, manager.bl)
    manager.align.assert_called_once_with(manager.gb)


def test_configure_preserves_disabled_stages(tmp_path: Path) -> None:
    config_file = tmp_path / "pipeline.yml"
    config_file.write_text("{}\n", encoding="utf-8")
    manager = DataMana()

    manager.configure(config_file)

    assert manager.pm is None
    assert manager.db is None
    assert manager.bl is None
    assert manager.gb is None
    assert manager.al is None


def test_database_dispatches_nested_and_flat_strategies() -> None:
    manager = DataMana()
    nested_action = mock.Mock()
    shared_action = mock.Mock()
    flat_action = mock.Mock()
    database_manager = SimpleNamespace(
        database_dict={
            "nested": [
                {"specific": nested_action, "shared": shared_action},
                {"specific": {"value": 1}, "shared": {"value": 2}},
            ],
            "flat": [flat_action, {"value": 3}],
        }
    )

    with mock.patch(
        "OrthoEvol.Manager.data_management.BaseDatabaseManagement",
        return_value=database_manager,
    ):
        manager.database(mock.sentinel.project, {"driver": "sqlite"})

    nested_action.assert_called_once_with(value=1)
    shared_action.assert_called_once_with(value=2)
    flat_action.assert_called_once_with(value=3)


def test_blast_builds_and_starts_configured_workflow() -> None:
    manager = DataMana()
    manager.Management_config = {"project": "example"}
    blast = mock.Mock()
    blast.blast_human = ["NM_000000"]

    with mock.patch(
        "OrthoEvol.Manager.data_management.OrthoBlastN",
        return_value=blast,
    ) as blast_class:
        manager.blast(mock.sentinel.project, {"method": 1})

    blast_class.assert_called_once_with(
        proj_mana=mock.sentinel.project,
        project="example",
        method=1,
    )
    blast.blast_config.assert_called_once_with(
        ["NM_000000"],
        "Homo_sapiens",
        auto_start=True,
    )


def test_genbank_uses_completed_blast_results() -> None:
    manager = DataMana()
    manager.Management_config = {"project": "example"}
    manager.GenBank_config = {"solo": True}
    blast = object.__new__(OrthoBlastN)
    blast.org_list = ["Homo_sapiens"]
    blast.gene_dict = {"GENE1": {"Homo_sapiens": "NM_000000"}}
    genbank = mock.Mock()

    with mock.patch(
        "OrthoEvol.Manager.data_management.GenBank",
        return_value=genbank,
    ):
        manager.genbank(mock.sentinel.project, blast)

    genbank.create_post_blast_gbk_records.assert_called_once_with(
        blast.org_list,
        blast.gene_dict,
    )


def test_genbank_fetches_configured_accessions_without_blast() -> None:
    manager = DataMana()
    manager.Management_config = {"project": "example"}
    manager.GenBank_config = {"solo": True}
    manager.CompGenAnalysis_config = {"taxon_file": "taxa.csv"}
    genbank = mock.Mock()
    comparative_analysis = SimpleNamespace(
        tier_frame_dict={"Tier1": SimpleNamespace(T=["GENE1"])},
        org_list=["Homo_sapiens"],
        gene_dict={"GENE1": {"Homo_sapiens": "nm_000000.2"}},
    )
    project_manager = SimpleNamespace(project="example")

    with mock.patch(
        "OrthoEvol.Manager.data_management.GenBank",
        return_value=genbank,
    ):
        with mock.patch(
            "OrthoEvol.Manager.data_management.BaseComparativeGenetics",
            return_value=comparative_analysis,
        ):
            manager.genbank(project_manager, None)

    genbank.get_gbk_file.assert_called_once_with(
        "NM_000000",
        "GENE1",
        "Homo_sapiens",
        server_flag=False,
    )


def test_align_dispatches_each_configured_aligner() -> None:
    manager = DataMana()
    manager.Management_config = {"project": "example"}
    manager.Alignment_config = {"Clustalo_config": {}}
    first_aligner = mock.Mock()
    second_aligner = mock.Mock()
    alignment = SimpleNamespace(
        alignment_dict={
            "first": [first_aligner, {"threads": 2}],
            "second": [second_aligner, {"iterations": 3}],
        }
    )

    with mock.patch(
        "OrthoEvol.Manager.data_management.MSA",
        return_value=alignment,
    ):
        manager.align(mock.sentinel.genbank)

    first_aligner.assert_called_once_with(threads=2)
    second_aligner.assert_called_once_with(iterations=3)
