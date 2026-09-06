"""Tests for the post-BLAST Excel report."""

from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from OrthoEvol.Orthologs.Blast.comparative_genetics import ComparativeGenetics


def build_comparative_genetics(tmp_path: Path) -> ComparativeGenetics:
    """Build only the state required to create a post-BLAST report."""
    analysis = object.__new__(ComparativeGenetics)
    analysis.data = tmp_path
    analysis.project = "example"
    analysis.postblastlog = Mock()
    analysis.dup_acc_count = {"NM_001": 2}
    analysis.dup_gene_count = {"GENE1": 1}
    analysis.duplicated_genes = {
        "GENE1": {"NM_001": ["Human", "Mouse"]}
    }
    analysis.dup_org_count = {"Human": 1}
    analysis.duplicated_organisms = {
        "Human": {"NM_001": ["GENE1", "GENE2"]}
    }
    analysis.duplicated_random = {
        "NM_002": [["GENE2", "Mouse"], ["GENE3", "Rat"]]
    }
    analysis.duplicated_other = {}
    analysis.missing_organsims = {
        "Mouse": {"missing genes": ["GENE3"], "count": 1}
    }
    analysis.missing_genes = {
        "GENE4": {"missing organisms": ["Rat"], "count": 1}
    }
    return analysis


def test_post_blast_analysis_writes_expected_workbook(tmp_path: Path) -> None:
    analysis = build_comparative_genetics(tmp_path)

    output_path = analysis.post_blast_analysis(["GENE5"])

    assert output_path == tmp_path / "example_postblastanalysis.xlsx"
    assert output_path.is_file()
    worksheets = pd.read_excel(output_path, sheet_name=None, index_col=0)
    assert list(worksheets) == [
        "Removed Genes",
        "Duplicate Count by Accession",
        "Duplicate Count by Gene",
        "Duplicate Org Groups by Gene",
        "Duplicate Count by Org",
        "Duplicate Gene Groups by Org",
        "Random Duplicates",
        "Missing Genes Count",
        "Missing Genes by Org",
        "Missing Organisms Count",
        "Missing Organisms by Gene",
    ]
    assert worksheets["Removed Genes"]["Removed Genes"].tolist() == ["GENE5"]
    assert worksheets["Duplicate Count by Accession"].at["NM_001", "Count"] == 2
    assert worksheets["Missing Genes by Org"].at["Mouse", 0] == "GENE3"
    assert worksheets["Missing Organisms by Gene"].at["GENE4", 0] == "Rat"
    analysis.postblastlog.info.assert_called_once_with(
        f"Post-BLAST analysis written to {output_path}."
    )


def test_post_blast_analysis_skips_empty_workbook(tmp_path: Path) -> None:
    analysis = build_comparative_genetics(tmp_path)
    analysis.dup_acc_count = {}
    analysis.dup_gene_count = {}
    analysis.duplicated_genes = {}
    analysis.dup_org_count = {}
    analysis.duplicated_organisms = {}
    analysis.duplicated_random = {}
    analysis.duplicated_other = {}
    analysis.missing_organsims = {}
    analysis.missing_genes = {}

    output_path = analysis.post_blast_analysis([])

    assert output_path is None
    assert not (tmp_path / "example_postblastanalysis.xlsx").exists()
    analysis.postblastlog.warning.assert_called_once_with(
        "Post-BLAST analysis contained no reportable results."
    )


def test_post_blast_analysis_rejects_malformed_missing_data(
    tmp_path: Path,
) -> None:
    analysis = build_comparative_genetics(tmp_path)
    analysis.missing_genes = {"GENE4": {"count": 1}}

    with pytest.raises(
        ValueError,
        match="GENE4.*missing organisms",
    ):
        analysis.post_blast_analysis([])
