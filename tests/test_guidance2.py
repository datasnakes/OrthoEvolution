"""Unit tests for GUIDANCE2 configuration and control flow."""

from pathlib import Path
from unittest import mock

import pytest

from OrthoEvol.Orthologs.Align.guidance2 import Guidance2Commandline
from OrthoEvol.Orthologs.Align.msa import MultipleSequenceAlignment


def build_alignment_without_external_tools(raw_data: Path) -> MultipleSequenceAlignment:
    """Construct the class around the pure helpers without loading project data."""
    alignment = object.__new__(MultipleSequenceAlignment)
    alignment.raw_data = raw_data
    return alignment


def test_guidance_paths_cover_amino_acid_column_filtering(
    tmp_path: Path,
) -> None:
    alignment = build_alignment_without_external_tools(tmp_path)

    paths = alignment._build_guidance_paths("gene.faa", "aa")

    assert paths.sequence_column_filtered == tmp_path / "gene/gene_G2sfcf_aa.aln"


def test_guidance_paths_reject_unsupported_sequence_type(tmp_path: Path) -> None:
    alignment = build_alignment_without_external_tools(tmp_path)

    with pytest.raises(ValueError, match="Unsupported GUIDANCE2 sequence type"):
        alignment._build_guidance_paths("gene.fna", "codon")


@pytest.mark.parametrize(
    ("removed_sequence_count", "iteration", "maximum_iterations", "expected"),
    [
        (0, 1, 5, True),
        (2, 1, 5, False),
        (2, 5, 5, True),
    ],
)
def test_guidance_stop_condition(
    removed_sequence_count: int,
    iteration: int,
    maximum_iterations: int,
    expected: bool,
) -> None:
    assert (
        MultipleSequenceAlignment._should_stop_guidance(
            removed_sequence_count,
            iteration,
            maximum_iterations,
        )
        is expected
    )


def test_sequence_filter_output_directory_is_stable(tmp_path: Path) -> None:
    output_directory = (
        MultipleSequenceAlignment._sequence_filter_output_directory(
            tmp_path / "gene",
            column_filter=None,
            mask_filter=None,
        )
    )

    assert output_directory == tmp_path / "gene/GUIDANCE2_sf"


def test_guidance_alignment_command_accepts_only_program_options(
    tmp_path: Path,
) -> None:
    sequence_file = tmp_path / "gene.faa"
    sequence_file.write_text(">gene\nM\n", encoding="utf-8")

    command = Guidance2Commandline(
        seqFile=sequence_file,
        msaProgram="MAFFT",
        seqType="aa",
        outDir=tmp_path / "output",
    )

    assert str(command).startswith("guidance --seqFile")


def test_guidance_mask_command_does_not_change_directory(tmp_path: Path) -> None:
    original_directory = Path.cwd()

    command = Guidance2Commandline(
        align=False,
        maskFile=tmp_path / "alignment.fasta",
        rprScores=tmp_path / "scores.txt",
        output=tmp_path / "masked.fasta",
        maskCutoff=0.6,
        seqType="aa",
    )

    assert str(command).startswith("maskLowScoreResidues")
    assert Path.cwd() == original_directory


def test_guidance_does_not_forward_iteration_controls(tmp_path: Path) -> None:
    alignment = build_alignment_without_external_tools(tmp_path)
    alignment.guidancelog = mock.Mock()
    sequence_file = tmp_path / "gene.faa"
    sequence_file.write_text(">gene\nM\n", encoding="utf-8")

    with mock.patch.object(alignment, "_run_guidance_command") as run_command:
        with mock.patch("OrthoEvol.Orthologs.Align.msa.shutil.copy"):
            alignment.guidance2(
                sequence_file,
                "MAFFT",
                "aa",
                columnFilter=0.9,
                iterations=3,
                increment=0.1,
            )

    command_options = run_command.call_args.kwargs
    assert "iterations" not in command_options
    assert "increment" not in command_options
