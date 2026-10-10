"""Focused tests for BLAST workflow control."""

import io
import sqlite3
import tarfile
from pathlib import Path
from subprocess import CalledProcessError
from types import SimpleNamespace
from unittest import mock
from xml.etree.ElementTree import ParseError

import pytest

from OrthoEvol.Orthologs.Blast.blast import ApplicationError, BaseBlastN, BlastFailure, OrthoBlastN
from OrthoEvol.Orthologs.Blast.comparative_genetics import BaseComparativeGenetics
from OrthoEvol.Tools.taxonomy import (
    prepare_taxonomy_database,
    resolve_species,
    taxonomy_database,
)


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


@pytest.mark.parametrize("method", [1, 2, None, ""])
def test_select_method_returns_query_configuration(method: int | str | None) -> None:
    parameters, query = BaseBlastN.select_method(method)

    assert parameters["query"] == ""
    assert query["temp fasta"] == ""
    assert parameters.get("remote") == ("True" if method == 2 else None)


def test_select_method_rejects_unknown_method() -> None:
    with pytest.raises(ValueError, match="not a blast method"):
        BaseBlastN.select_method(3)


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


def test_create_temp_fasta_propagates_failure_with_spaced_path(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    query_config = {
        "query": "NM_001",
        "db": "refseq",
        "temp fasta": str(tmp_path / "query with spaces.fasta"),
    }

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.run",
        side_effect=CalledProcessError(1, "blastdbcmd", stderr=b"missing query"),
    ) as command, pytest.raises(CalledProcessError):
        blast._create_temp_fasta("NM_001", "GeneA", query_config)

    assert command.call_args.args[0][-1] == query_config["temp fasta"]
    assert command.call_args.kwargs["check"] is True


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


@pytest.mark.parametrize(
    "error",
    [ParseError("invalid XML"), ApplicationError(1, "blastn", stderr="bad database")],
)
def test_blastn_wrapper_removes_failed_output_and_raises(
    tmp_path: Path, error: Exception
) -> None:
    blast = build_blast(tmp_path)
    blast.get_time = mock.Mock(return_value=1.0)
    command = mock.Mock(side_effect=error)
    xml_path = tmp_path / "result.xml"

    with mock.patch(
        "OrthoEvol.Orthologs.Blast.blast.NcbiblastnCommandline",
        return_value=command,
    ), pytest.raises(type(error)):
        blast.blastn_wrapper(
            "GeneA", "Mus_musculus", {}, str(xml_path), tmp_path
        )

    assert not xml_path.exists()


def test_external_accession_path_keeps_derived_files_in_project(tmp_path: Path) -> None:
    source = tmp_path / "external.v1.csv"
    source.write_text("Tier,Gene,Homo_sapiens\n1,ADRA1A,NM_000680.3\n")
    analysis = object.__new__(BaseComparativeGenetics)
    analysis.acc_file = source
    analysis.project = "analysis"
    analysis.project_index = tmp_path / "index"
    analysis.project_index.mkdir()
    analysis.data = tmp_path / "data"

    analysis._copy_accession_file(copy_from_package=False)
    analysis._initialize_accession_paths(source)

    assert analysis.acc_csv_path.read_text() == source.read_text()
    assert analysis.acc_csv_path.parent == analysis.project_index
    assert analysis.building_file_path == analysis.data / "external.v1building.csv"


def test_runblast_does_not_complete_after_invalid_cached_xml(tmp_path: Path) -> None:
    blast = build_blast(tmp_path)
    blast.date_format = "test date"
    blast.org_list = ["Homo_sapiens", "Mus_musculus"]
    blast.taxon_dict = {"Mus_musculus": 10090}
    blast.create_maf = mock.Mock()
    xml_path = tmp_path / "GeneA" / "BLAST" / "GeneA_Mus_musculus.xml"
    xml_path.parent.mkdir(parents=True)
    xml_path.write_text("<broken>")

    with pytest.raises(ParseError):
        blast.runblast(
            genes=["GeneA"], query_organism="Homo_sapiens", pre_configured=True
        )

    assert not xml_path.exists()
    blast.create_maf.assert_not_called()


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


@pytest.fixture
def taxonomy_path(tmp_path: Path) -> Path:
    """Build a real tiny NCBI dump including the duplicate synonym that broke ETE."""
    species = [(1, 1, "root"), (9606, 1, "Homo sapiens"), (10090, 1, "Mus musculus")]
    names = "".join(
        f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
        for taxid, _, name in species
    )
    names += "9606\t|\tHuman\t|\t\t|\tsynonym\t|\n9606\t|\thuman\t|\t\t|\tsynonym\t|\n"
    nodes = "".join(
        f"{taxid}\t|\t{parent}\t|\t{'no rank' if taxid == 1 else 'species'}\t|\n"
        for taxid, parent, _ in species
    )
    archive = tmp_path / "taxdump.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        for name, content in [
            ("names.dmp", names),
            ("nodes.dmp", nodes),
            ("merged.dmp", "99999\t|\t9606\t|\n"),
        ]:
            member = tarfile.TarInfo(name)
            member.size = len(content.encode())
            bundle.addfile(member, io.BytesIO(content.encode()))
    return prepare_taxonomy_database(
        tmp_path / "taxa.sqlite",
        taxdump_file=archive,
        required_species=["Homo_sapiens", "Mus_musculus"],
    )


def test_taxonomy_failed_replacement_preserves_existing_database(
    taxonomy_path: Path,
) -> None:
    before = taxonomy_path.read_bytes()
    cache = Path(f"{taxonomy_path}.traverse.pkl").read_bytes()
    with pytest.raises(ValueError, match="Missing or ambiguous"):
        prepare_taxonomy_database(
            taxonomy_path,
            taxdump_file=taxonomy_path.parent / "taxdump.tar.gz",
            required_species=["Unknown_species"],
            overwrite=True,
        )
    assert taxonomy_path.read_bytes() == before
    assert Path(f"{taxonomy_path}.traverse.pkl").read_bytes() == cache
    with taxonomy_database(taxonomy_path) as taxonomy:
        assert resolve_species(taxonomy, ["human"]) == {"human": 9606}


@pytest.mark.parametrize("state", ["missing", "empty", "partial"])
def test_invalid_taxonomy_never_triggers_download(tmp_path: Path, state: str) -> None:
    path = tmp_path / "invalid.sqlite"
    if state != "missing":
        with sqlite3.connect(path) as db:
            db.executescript(
                "CREATE TABLE stats(version); INSERT INTO stats VALUES(2); CREATE TABLE species(taxid, track);"
            )
            if state == "partial":
                db.execute("INSERT INTO species VALUES(1, NULL)")
    with mock.patch("ete3.NCBITaxa.update_taxonomy_database") as update:
        with pytest.raises((FileNotFoundError, ValueError)):
            with taxonomy_database(path):
                pass
        update.assert_not_called()


@pytest.mark.parametrize("row", ["1,ADRA1A,,NM_001\n1,ADRA1A,,NM_002"])
def test_invalid_reference_table_fails_before_taxonomy(
    tmp_path: Path, row: str
) -> None:
    source = tmp_path / "accessions.csv"
    source.write_text(f"Tier,Gene,Mus_musculus,Homo_sapiens\n{row}\n")
    with pytest.raises(ValueError, match="Missing query accessions|duplicate genes"):
        OrthoBlastN(
            project="test",
            project_path=tmp_path,
            acc_file=source,
            copy_from_package=False,
            taxonomy_db=tmp_path / "missing.sqlite",
        )
