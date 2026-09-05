"""Focused tests for GenBank record retrieval and validation."""

from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from OrthoEvol.Orthologs.GenBank.genbank import GenBank


def build_genbank(tmp_path: Path) -> GenBank:
    """Build only the state used by the tested GenBank methods."""
    genbank = object.__new__(GenBank)
    genbank.raw_data = tmp_path / "raw_data"
    genbank.ncbi_db_repo = tmp_path / "databases"
    genbank.db_files_list = ["records.db"]
    genbank.duplicated_dict = {}
    genbank.genbanklog = mock.Mock()
    return genbank


def test_create_post_blast_records_normalizes_accessions(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    genbank.tier_frame_dict = {1: SimpleNamespace(T=["gene_a"])}
    genbank.get_gbk_file = mock.Mock()

    genbank.create_post_blast_gbk_records(
        ["Mus_musculus"],
        {"gene_a": {"Mus_musculus": "nm_001.4"}},
    )

    genbank.get_gbk_file.assert_called_once_with(
        "NM_001", "gene_a", "Mus_musculus", server_flag=False
    )


def test_get_gbk_file_writes_first_matching_record(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    record = mock.Mock()
    record.format.return_value = "LOCUS record\n"
    database = mock.Mock()
    database.lookup.return_value = record
    server = mock.MagicMock()
    server.keys.return_value = ["refseq"]
    server.__getitem__.return_value = database
    genbank.gbk_quality_control = mock.Mock()

    with mock.patch(
        "OrthoEvol.Orthologs.GenBank.genbank.BioSeqDatabase.open_database",
        return_value=server,
    ):
        genbank.get_gbk_file("NM_001", "GeneA", "Mus_musculus")

    output = tmp_path / "raw_data" / "GeneA" / "GENBANK" / "GeneA_Mus_musculus.gbk"
    assert output.read_text(encoding="utf-8") == "LOCUS record\n"
    genbank.gbk_quality_control.assert_called_once_with(
        output, "GeneA", "Mus_musculus"
    )


def test_get_gbk_file_raises_when_accession_is_absent(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    database = mock.Mock()
    database.lookup.side_effect = IndexError
    server = mock.MagicMock()
    server.keys.return_value = ["refseq"]
    server.__getitem__.return_value = database

    with mock.patch(
        "OrthoEvol.Orthologs.GenBank.genbank.BioSeqDatabase.open_database",
        return_value=server,
    ):
        with pytest.raises(FileNotFoundError):
            genbank.get_gbk_file("MISSING", "GeneA", "Mus_musculus")


@pytest.mark.parametrize(
    ("record_gene", "record_organism", "raises"),
    [
        ("GeneA", "Mus musculus", False),
        ("GeneB", "Mus musculus", True),
        ("GeneA", "Rattus norvegicus", True),
    ],
)
def test_gbk_quality_control_validates_gene_and_organism(
    record_gene: str,
    record_organism: str,
    raises: bool,
    tmp_path: Path,
) -> None:
    genbank = build_genbank(tmp_path)
    record = SimpleNamespace(
        id="NM_001.4",
        features=[
            SimpleNamespace(qualifiers={"organism": [record_organism]}),
            SimpleNamespace(qualifiers={"gene": [record_gene]}),
        ],
    )

    with mock.patch(
        "OrthoEvol.Orthologs.GenBank.genbank.SeqIO.read", return_value=record
    ):
        if raises:
            with pytest.raises(BrokenPipeError):
                genbank.gbk_quality_control(
                    tmp_path / "record.gbk", "GeneA", "Mus_musculus"
                )
        else:
            genbank.gbk_quality_control(
                tmp_path / "record.gbk", "GeneA", "Mus_musculus"
            )

    if not raises:
        assert genbank.duplicated_dict["validated"] == {
            "NM_001.4": ["GeneA", "Mus_musculus"]
        }


def test_write_fasta_files_writes_coding_sequences(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    genbank.solo = True
    genbank.multi = True
    genbank.min_fasta = True
    feature = SimpleNamespace(
        type="CDS",
        qualifiers={
            "GI:456": [],
            "protein_id": ["NP_001"],
            "product": ["protein A"],
            "translation": ["MK"],
            "note": ["coding sequence"],
        },
        extract=mock.Mock(return_value="ATGAAG"),
    )
    record = SimpleNamespace(
        id="NM_001",
        description="GeneA transcript",
        annotations={"gi": "123"},
        seq="ATGAAG",
        features=[feature],
    )

    genbank.write_fasta_files(
        record, {"NM_001": ["GeneA", "Mus_musculus"]}
    )

    output_dir = tmp_path / "raw_data" / "GeneA" / "GENBANK"
    assert (output_dir / "GeneA_Mus_musculusCDS.ffn").read_text() == ">M_musculus\nATGAAG\n"
    assert (output_dir / "GeneACDS.faa").read_text() == ">M_musculus\nMK\n"


@pytest.mark.parametrize(
    ("feature_type", "expected_name"),
    [
        ("misc_feature", "GeneA_Mus_musculus_misc_feature.fna"),
        ("source", "GeneA_Mus_musculus_source.fasta"),
    ],
)
def test_fasta_writers_handle_non_coding_features(
    feature_type: str, expected_name: str, tmp_path: Path
) -> None:
    genbank = build_genbank(tmp_path)
    output_dir = tmp_path / "fasta"
    formatter = {
        "feat_type": feature_type,
        "feat_type_rank": feature_type,
        "path": str(output_dir),
        "gene": "GeneA",
        "org": "Mus_musculus",
        "na_gi": "123",
        "na_acc_n": "NM_001",
        "na_description": "GeneA transcript",
        "na_misc_feat": "conserved region",
        "na_seq": "ATGAAG",
    }

    genbank.solo_fasta(">sequence\nATGAAG\n", ">protein\nMK\n", formatter)
    if feature_type == "misc_feature":
        genbank.multi_fasta(">sequence\nATGAAG\n", ">protein\nMK\n", formatter)

    assert (output_dir / expected_name).is_file()


def test_get_fasta_files_reads_biosql_records(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    genbank.target_gbk_db_path = tmp_path / "project_databases"
    genbank.target_gbk_db_path.mkdir()
    (genbank.target_gbk_db_path / "records.db").touch()
    record = object()
    database = mock.MagicMock()
    database.keys.return_value = ["NM_001"]
    database.lookup.return_value = record
    server = mock.MagicMock()
    server.keys.return_value = ["GeneA"]
    server.__getitem__.return_value = database
    genbank.write_fasta_files = mock.Mock()

    with mock.patch(
        "OrthoEvol.Orthologs.GenBank.genbank.BioSeqDatabase.open_database",
        return_value=server,
    ):
        genbank.get_fasta_files({"NM_001": ["GeneA", "Mus_musculus"]})

    genbank.write_fasta_files.assert_called_once_with(
        record, {"NM_001": ["GeneA", "Mus_musculus"]}
    )


def test_get_fasta_files_reads_local_genbank_records(tmp_path: Path) -> None:
    genbank = build_genbank(tmp_path)
    genbank.target_gbk_files_path = tmp_path / "records"
    genbank.target_gbk_files_path.mkdir()
    (genbank.target_gbk_files_path / "record.gbk").touch()
    record = object()
    genbank.write_fasta_files = mock.Mock()

    with mock.patch(
        "OrthoEvol.Orthologs.GenBank.genbank.SeqIO.read", return_value=record
    ):
        genbank.get_fasta_files(
            {"NM_001": ["GeneA", "Mus_musculus"]}, db=False
        )

    genbank.write_fasta_files.assert_called_once_with(
        record, {"NM_001": ["GeneA", "Mus_musculus"]}
    )
