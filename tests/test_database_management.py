"""Tests for database strategy construction without external services."""

from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Manager.database_management import DatabaseManagement


def build_database_manager(tmp_path: Path) -> DatabaseManagement:
    """Build only the state required by the strategy helper methods."""
    manager = object.__new__(DatabaseManagement)
    manager.user_db = tmp_path / "databases"
    manager.user_archive = tmp_path / "archive"
    manager.db_mana_utils = SimpleNamespace(archive=Mock(name="archive"))
    manager.download_blast_database = Mock(name="download_blast_database")
    manager.download_windowmasker_files = Mock(name="download_windowmasker_files")
    manager.download_ncbi_taxonomy_dump_files = Mock(
        name="download_ncbi_taxonomy_dump_files"
    )
    manager.gene_data = SimpleNamespace(taxon_ids=[9606, 9544])
    return manager


def test_prepare_child_strategies_does_not_mutate_callers() -> None:
    original_strategy = {"custom": "value"}

    (prepared_strategy,) = DatabaseManagement._prepare_child_strategies(
        (original_strategy,),
        configure_flag=True,
        archive_flag=True,
        delete_flag=True,
    )

    assert original_strategy == {"custom": "value"}
    assert prepared_strategy == {
        "custom": "value",
        "configure_flag": True,
        "archive_flag": None,
        "delete_flag": None,
    }


def test_refseq_upload_reports_retired_scheduler(tmp_path: Path) -> None:
    manager = build_database_manager(tmp_path)

    with pytest.raises(
        OrthoEvolDeprecationWarning,
        match="retired SGE subsystem",
    ):
        DatabaseManagement.NCBI_refseq_release(manager, upload_flag=True)


def test_ncbi_preserves_dispatch_order_without_mutating_inputs(
    tmp_path: Path,
) -> None:
    manager = build_database_manager(tmp_path)
    manager.NCBI_blast = Mock(
        return_value=(
            OrderedDict({"NCBI_blast": ["blast_action"]}),
            OrderedDict({"NCBI_blast": [{}]}),
        )
    )
    manager.NCBI_pub_taxonomy = Mock(
        return_value=(
            OrderedDict({"NCBI_pub_taxonomy": ["taxonomy_action"]}),
            OrderedDict({"NCBI_pub_taxonomy": [{}]}),
        )
    )
    manager.NCBI_refseq_release = Mock(
        return_value=(
            OrderedDict({"NCBI_refseq_release": ["refseq_action"]}),
            OrderedDict({"NCBI_refseq_release": [{}]}),
        )
    )
    blast_options: dict[str, object] = {}
    taxonomy_options: dict[str, object] = {}
    refseq_options: dict[str, object] = {}

    dispatcher, configuration = DatabaseManagement.NCBI(
        manager,
        blast_options,
        taxonomy_options,
        refseq_options,
        configure_flag=True,
        archive_flag=True,
        delete_flag=True,
    )

    assert blast_options == {}
    assert taxonomy_options == {}
    assert refseq_options == {}
    assert list(dispatcher) == [
        "NCBI",
        "NCBI_blast",
        "NCBI_pub_taxonomy",
        "NCBI_refseq_release",
    ]
    assert list(configuration) == list(dispatcher)
    assert configuration["NCBI"] == [
        {
            "database_path": str(manager.user_db),
            "archive_path": str(manager.user_archive),
            "option": "NCBI",
            "delete_flag": True,
        }
    ]
    manager.NCBI_blast.assert_called_once_with(
        configure_flag=True,
        archive_flag=None,
        delete_flag=None,
    )


def test_strategy_dispatcher_maps_every_strategy_and_preserves_full_reset(
    tmp_path: Path,
) -> None:
    manager = build_database_manager(tmp_path)
    strategy_methods = {
        "NCBI": "NCBI",
        "NCBI_blast": "NCBI_blast",
        "Full": "full",
        "NCBI_blast_db": "ncbi_blast_db",
        "NCBI_blast_windowmaskerfiles": "ncbi_blast_windowmasker_files",
        "NCBI_pub_taxonomy": "NCBI_pub_taxonomy",
        "NCBI_refseq_release": "NCBI_refseq_release",
        "ITIS": "itis",
        "ITIS_taxonomy": "itis_taxonomy",
    }
    for strategy_name, method_name in strategy_methods.items():
        result = (
            OrderedDict({strategy_name: [f"{strategy_name}_action"]}),
            OrderedDict({strategy_name: [{"strategy": strategy_name}]}),
        )
        setattr(manager, method_name, Mock(return_value=result))

    configuration = OrderedDict(
        (strategy_name, {"value": strategy_name})
        for strategy_name in (*strategy_methods, "unknown")
    )
    dispatcher, strategy_config = manager.get_strategy_dispatcher(configuration)

    expected_strategies = list(strategy_methods)[2:]
    assert list(dispatcher) == expected_strategies
    assert list(strategy_config) == expected_strategies
    for strategy_name, method_name in strategy_methods.items():
        getattr(manager, method_name).assert_called_once_with(value=strategy_name)


@pytest.mark.parametrize(
    (
        "method_name",
        "strategy_name",
        "method_kwargs",
        "configure_action_name",
        "configure_kwargs",
        "archive_option",
        "expected_database_path",
    ),
    (
        (
            "ncbi_blast_db",
            "NCBI_blast_db",
            {},
            "download_blast_database",
            {"database_name": "refseq_rna"},
            "NCBI_blast",
            None,
        ),
        (
            "ncbi_blast_windowmasker_files",
            "NCBI_blast_windowmasker_files",
            {"taxonomy_ids": [1]},
            "download_windowmasker_files",
            {"taxonomy_ids": [9606, 9544]},
            "NCBI_blast_windowmasker_files",
            "default",
        ),
        (
            "NCBI_pub_taxonomy",
            "NCBI_pub_taxonomy",
            {},
            "download_ncbi_taxonomy_dump_files",
            {},
            "NCBI_pub_taxonomy",
            "default",
        ),
    ),
)
def test_leaf_strategies_preserve_action_order_and_configuration(
    tmp_path: Path,
    method_name: str,
    strategy_name: str,
    method_kwargs: dict[str, object],
    configure_action_name: str,
    configure_kwargs: dict[str, object],
    archive_option: str,
    expected_database_path: str | None,
) -> None:
    manager = build_database_manager(tmp_path)
    strategy = getattr(manager, method_name)

    dispatcher, configuration = strategy(
        configure_flag=True,
        archive_flag=True,
        delete_flag=True,
        **method_kwargs,
    )

    database_path = (
        str(manager.user_db)
        if expected_database_path == "default"
        else expected_database_path
    )
    assert list(dispatcher) == [strategy_name]
    assert dispatcher[strategy_name] == [
        manager.db_mana_utils.archive,
        getattr(manager, configure_action_name),
    ]
    assert configuration[strategy_name] == [
        {
            "database_path": database_path,
            "archive_path": str(manager.user_archive),
            "option": archive_option,
            "delete_flag": True,
        },
        configure_kwargs,
    ]
