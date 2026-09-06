"""Tests for duplicate accession classification and summaries."""

import pytest

from OrthoEvol.utilities import FullUtilities


def test_duplicate_analysis_classifies_each_relationship() -> None:
    accession_data = {
        "SAME_ORGANISM": [["G1", "O1"], ["G2", "O1"]],
        "SAME_ORGANISM_2": [["G2", "O1"], ["G3", "O1"]],
        "SAME_GENE": [["G1", "O1"], ["G1", "O2"]],
        "SAME_GENE_2": [["G1", "O2"], ["G1", "O3"]],
        "RANDOM": [["G1", "O1"], ["G2", "O2"]],
        "MIXED": [["G1", "O1"], ["G1", "O2"], ["G2", "O2"]],
        "UNIQUE": [["G2", "O1"]],
    }

    analysis = FullUtilities().analyze_duplicate_accessions(
        accession_data,
        gene_list=["G1", "G2", "G3"],
        org_list=["O1", "O2", "O3"],
    )

    assert set(analysis.groups["accessions"]) == {
        "SAME_ORGANISM",
        "SAME_ORGANISM_2",
        "SAME_GENE",
        "SAME_GENE_2",
        "RANDOM",
        "MIXED",
    }
    assert analysis.groups["organisms"] == {
        "O1": {
            "SAME_ORGANISM": ["G1", "G2"],
            "SAME_ORGANISM_2": ["G2", "G3"],
        }
    }
    assert analysis.groups["genes"] == {
        "G1": {
            "SAME_GENE": ["O1", "O2"],
            "SAME_GENE_2": ["O2", "O3"],
        }
    }
    assert analysis.groups["random"] == {
        "RANDOM": [["G1", "O1"], ["G2", "O2"]]
    }
    assert analysis.groups["other"] == {
        "MIXED": [["G1", "O1"], ["G1", "O2"], ["G2", "O2"]]
    }
    assert analysis.accession_counts == {
        "SAME_ORGANISM": 2,
        "SAME_ORGANISM_2": 2,
        "SAME_GENE": 2,
        "SAME_GENE_2": 2,
        "RANDOM": 2,
        "MIXED": 3,
    }
    assert analysis.gene_counts == {"G1": 2, "G2": 0, "G3": 0}
    assert analysis.organism_counts == {"O1": 2, "O2": 0, "O3": 0}


def test_get_dup_acc_preserves_dictionary_contract() -> None:
    duplicate_groups = FullUtilities().get_dup_acc(
        {"A1": [["G1", "O1"], ["G1", "O2"]]},
        gene_list=["G1"],
        org_list=["O1", "O2"],
    )

    assert isinstance(duplicate_groups, dict)
    assert list(duplicate_groups) == [
        "accessions",
        "genes",
        "organisms",
        "random",
        "other",
    ]


def test_duplicate_analysis_rejects_malformed_pairs() -> None:
    with pytest.raises(
        ValueError,
        match="BROKEN.*gene-organism pairs",
    ):
        FullUtilities().analyze_duplicate_accessions(
            {"BROKEN": [["G1"]]},
            gene_list=["G1"],
            org_list=["O1"],
        )
