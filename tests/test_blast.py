"""Focused tests for BLAST workflow control."""

from pathlib import Path
from types import SimpleNamespace
from unittest import mock
from xml.etree.ElementTree import ParseError

import pytest

from OrthoEvol.Orthologs.Blast.blast import BaseBlastN, BlastFailure


def build_blast(tmp_path: Path, method: int | None = 1) -> BaseBlastN:
    """Build only the state required by workflow methods."""
    blast = object.__new__(BaseBlastN)
    blast.method = method
    blast.raw_data = tmp_path
    blast.blastn_log = mock.Mock()
    blast.blast_utils = mock.Mock()
    blast.blastn_parameters = {}
    blast.query_config = {"query": "", "db": "refseq", "temp fasta": ""}
    blast.removed_genes = []
    return blast


@pytest.mark.parametrize("method", [1, 2, None])
def test_select_method_returns_query_configuration(method: int | None) -> None:
    parameters, query = BaseBlastN.select_method(method)

    assert parameters["query"] == ""
    assert query["temp fasta"] == ""
    assert parameters["db"] == "refseq_rna"
    assert parameters.get("remote") is (True if method == 2 else None)


@pytest.mark.parametrize("method", [3, ""])
def test_select_method_rejects_unknown_method(method: int | str) -> None:
    with pytest.raises(ValueError, match="not a blast method"):
        BaseBlastN.select_method(method)  # type: ignore[arg-type]


@pytest.mark.parametrize("method", [1, None])
def test_local_blast_requires_database_environment(method: int | None) -> None:
    with pytest.raises(EnvironmentError, match="BLASTDB is required"):
        BaseBlastN._validate_database_environment(method, {})


def test_remote_blast_does_not_require_database_environment() -> None:
    BaseBlastN._validate_database_environment(2, {})


def test_make_blast_dir_creates_parents(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    output = tmp_path / "GeneA" / "BLAST"

    blast._make_blast_dir("GeneA", output)

    assert output.is_dir()


def test_configure_resumes_remaining_queries(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.ref_species = None
    blast.building_file_path = tmp_path / "building.csv"
    blast.data = tmp_path
    blast.gene_list = ["GeneA", "GeneB"]
    blast.blast_human = ["NM_A", "NM_B"]
    blast.taxon_dict = {}
    blast.acc_dict = {"NM_A": [["GeneA"]], "NM_B": [["GeneB"]]}
    blast.blast_utils.gene_list_config.return_value = ["GeneB"]
    blast._make_blast_dir = mock.Mock()
    blast._create_temp_fasta = mock.Mock()

    blast.configure(blast.blast_human, "Homo_sapiens")

    assert blast.current_gene_list == ["GeneB"]
    blast._create_temp_fasta.assert_called_once_with(
        query="NM_B",
        gene="GeneB",
        query_config={
            "query": "NM_B",
            "db": "refseq",
            "temp fasta": str(tmp_path / "GeneB" / "BLAST" / "temp.fasta"),
        },
    )


def test_create_maf_copies_completed_files(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.building_file_path = tmp_path / "building.csv"
    blast.building_time_file_path = tmp_path / "building_time.csv"
    blast.complete_file_path = tmp_path / "complete.csv"
    blast.complete_time_file_path = tmp_path / "complete_time.csv"
    blast.building_file_path.write_text("accessions", encoding="utf-8")
    blast.building_time_file_path.write_text("timing", encoding="utf-8")

    blast.create_maf()

    assert blast.complete_file_path.read_text(encoding="utf-8") == "accessions"
    assert blast.complete_time_file_path.read_text(encoding="utf-8") == "timing"


def test_create_maf_reports_missing_building_files(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.building_file_path = tmp_path / "missing.csv"
    blast.building_time_file_path = tmp_path / "missing_time.csv"
    blast.complete_file_path = tmp_path / "complete.csv"
    blast.complete_time_file_path = tmp_path / "complete_time.csv"

    with pytest.raises(BlastFailure):
        blast.create_maf()


def test_parse_xml_selects_highest_matching_hit(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.add_accession = mock.Mock()
    blast.blast_utils.map_func = mock.Mock(side_effect=lambda hit: hit)
    low_hsp = SimpleNamespace(bitscore_raw=10)
    high_hsp = SimpleNamespace(bitscore_raw=20)
    hits = [
        SimpleNamespace(
            id="NM_LOW", id1="NM_LOW", id2="1", description="GeneA", hsps=[low_hsp]
        ),
        SimpleNamespace(
            id="NM_HIGH", id1="NM_HIGH", id2="2", description="GeneA", hsps=[high_hsp]
        ),
    ]
    query_result = mock.Mock()
    query_result.hit_map.return_value = hits
    xml_path = tmp_path / "result.xml"
    xml_path.write_text("<xml />", encoding="utf-8")

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.SearchIO.read", return_value=query_result
    ):
        blast.parse_xml(str(xml_path), "GeneA", "Mus_musculus")

    blast.add_accession.assert_called_once_with("GeneA", "Mus_musculus", "NM_HIGH")


@pytest.mark.parametrize("method", [1, 2])
def test_blastn_wrapper_writes_and_parses_output(
    method: int, tmp_path: Path
) -> None:
    blast = build_blast(tmp_path, method=method)
    blast.get_time = mock.Mock(side_effect=[1.0, 3.0])
    blast.add_blast_time = mock.Mock()
    blast.parse_xml = mock.Mock()
    command = mock.Mock(return_value=("<xml />", ""))
    xml_path = tmp_path / "result.xml"

    with (
        mock.patch(
            "OrthoEvol.Orthologs.Blast.blast.NcbiblastnCommandline",
            return_value=command,
        ),
        mock.patch("OrthoEvol.Orthologs.Blast.blast.time.sleep") as sleep,
    ):
        blast.blastn_wrapper(
            "GeneA",
            "Mus_musculus",
            {"query": "temp.fasta"},
            str(xml_path),
            tmp_path,
        )

    assert xml_path.read_text(encoding="utf-8") == "<xml />"
    blast.add_blast_time.assert_called_once_with(
        "GeneA", "Mus_musculus", 1.0, 3.0
    )
    blast.parse_xml.assert_called_once_with(
        str(xml_path), "GeneA", "Mus_musculus"
    )
    assert sleep.called is (method == 2)


@pytest.mark.parametrize("return_code", [0, 1])
def test_create_temp_fasta_handles_command_status(
    return_code: int, tmp_path: Path
) -> None:
    blast = build_blast(tmp_path)
    blast.gene_list = ["GeneA"]
    blast.current_gene_list = ["GeneA"]
    query_config = {
        "query": "NM_001",
        "db": "refseq",
        "temp fasta": str(tmp_path / "temp.fasta"),
    }

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.run",
        return_value=SimpleNamespace(returncode=return_code),
    ):
        blast._create_temp_fasta("NM_001", "GeneA", query_config)

    if return_code == 0:
        assert blast.gene_list == ["GeneA"]
        assert blast.removed_genes == []
    else:
        assert blast.gene_list == []
        assert blast.current_gene_list == []
        assert blast.removed_genes == ["GeneA"]


def test_runblast_dispatches_one_local_search(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.date_format = "test date"
    blast.org_list = ["Homo_sapiens", "Mus_musculus"]
    blast.taxon_dict = {"Mus_musculus": 10090}
    blast.tier_dict = {"GeneA": 1}
    blast.save_data = False
    blast.blastn_wrapper = mock.Mock()
    blast.create_maf = mock.Mock()
    gene_path = tmp_path / "GeneA" / "BLAST"
    gene_path.mkdir(parents=True)

    blast.runblast(
        genes=["GeneA"],
        query_organism="Homo_sapiens",
        pre_configured=True,
    )

    blast.blastn_wrapper.assert_called_once_with(
        xml_path=str(gene_path / "GeneA_Mus_musculus.xml"),
        gene_path=gene_path,
        parameters={
            "query": str(gene_path / "temp.fasta"),
            "taxids": 10090,
        },
        gene="GeneA",
        organism="Mus_musculus",
    )
    blast.create_maf.assert_called_once_with()


def test_blastn_wrapper_removes_unparseable_output(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.get_time = mock.Mock(return_value=1.0)
    command = mock.Mock(side_effect=ParseError("invalid XML"))
    xml_path = tmp_path / "result.xml"

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.NcbiblastnCommandline",
        return_value=command,
    ):
        blast.blastn_wrapper(
            "GeneA", "Mus_musculus", {}, str(xml_path), tmp_path
        )

    assert not xml_path.exists()


def test_parse_xml_skips_predicted_hit_and_lowercases_nonmatching_hit(
    tmp_path: Path,
) -> None:
    blast = build_blast(tmp_path)
    blast.add_accession = mock.Mock()
    blast.blast_utils.map_func = mock.Mock(side_effect=lambda hit: hit)
    hits = [
        SimpleNamespace(
            id="XR_001",
            id1="XR_001",
            id2="1",
            description="GeneA predicted RNA",
            hsps=[SimpleNamespace(bitscore_raw=30)],
        ),
        SimpleNamespace(
            id="NM_002",
            id1="NM_002",
            id2="2",
            description="unrelated transcript",
            hsps=[SimpleNamespace(bitscore_raw=20)],
        ),
    ]
    query_result = mock.Mock()
    query_result.hit_map.return_value = hits
    xml_path = tmp_path / "result.xml"
    xml_path.write_text("<xml />", encoding="utf-8")

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.SearchIO.read", return_value=query_result
    ):
        blast.parse_xml(str(xml_path), "GeneA", "Mus_musculus")

    blast.add_accession.assert_called_once_with(
        "GeneA", "Mus_musculus", "nm_002"
    )
