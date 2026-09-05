"""Tests for the Slurm command adapter."""

import subprocess
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from OrthoEvol.Tools.slurm import (SlurmClient, SlurmCommandNotFoundError,
                                   SlurmJob)


@patch("OrthoEvol.Tools.slurm.client.shutil.which", return_value="/usr/bin/sbatch")
@patch("OrthoEvol.Tools.slurm.client.subprocess.run")
def test_submit_returns_parsable_job_id(
    mock_run: Mock,
    _mock_which: Mock,
    tmp_path: Path,
) -> None:
    script = tmp_path / "analysis.sh"
    script.write_text("#!/bin/bash\ntrue\n", encoding="utf-8")
    mock_run.return_value = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout="12345;cheaha\n",
        stderr="",
    )

    job_id = SlurmClient().submit(script)

    assert job_id == "12345"
    mock_run.assert_called_once_with(
        ["sbatch", "--parsable", str(script)],
        check=True,
        capture_output=True,
        text=True,
    )


@patch("OrthoEvol.Tools.slurm.client.shutil.which", return_value="/usr/bin/squeue")
@patch("OrthoEvol.Tools.slurm.client.subprocess.run")
def test_active_jobs_parses_explicit_squeue_fields(
    mock_run: Mock,
    _mock_which: Mock,
) -> None:
    mock_run.return_value = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout=(
            "12345 | orthologs | RUNNING | medium | 00:03:21 | 1 | c0123\n"
            "12346 | alignment | PENDING | medium | 0:00 | 2 | Resources\n"
        ),
        stderr="",
    )

    jobs = SlurmClient().active_jobs(user="researcher")

    assert jobs == [
        SlurmJob(
            job_id="12345",
            name="orthologs",
            state="RUNNING",
            partition="medium",
            elapsed="00:03:21",
            nodes=1,
            node_list_or_reason="c0123",
        ),
        SlurmJob(
            job_id="12346",
            name="alignment",
            state="PENDING",
            partition="medium",
            elapsed="0:00",
            nodes=2,
            node_list_or_reason="Resources",
        ),
    ]
    mock_run.assert_called_once_with(
        [
            "squeue",
            "--noheader",
            "--user=researcher",
            "--format=%i|%j|%T|%P|%M|%D|%R",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


@patch("OrthoEvol.Tools.slurm.client.shutil.which", return_value="/usr/bin/sacct")
@patch("OrthoEvol.Tools.slurm.client.subprocess.run")
def test_job_history_parses_sacct_parsable_output(
    mock_run: Mock,
    _mock_which: Mock,
) -> None:
    mock_run.return_value = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout="12345|orthologs|COMPLETED|medium|00:04:12|1|c0123\n",
        stderr="",
    )

    jobs = SlurmClient().job_history("12345")

    assert jobs[0].state == "COMPLETED"
    assert jobs[0].nodes == 1
    mock_run.assert_called_once_with(
        [
            "sacct",
            "--allocations",
            "--noheader",
            "--parsable2",
            "--jobs=12345",
            "--format=JobIDRaw,JobName,State,Partition,Elapsed,AllocNodes,NodeList",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def test_submit_rejects_missing_script(tmp_path: Path) -> None:
    missing_script = tmp_path / "missing.sh"

    with pytest.raises(FileNotFoundError, match="missing.sh"):
        SlurmClient().submit(missing_script)


@patch("OrthoEvol.Tools.slurm.client.shutil.which", return_value=None)
def test_reports_missing_slurm_command(_mock_which: Mock) -> None:
    with pytest.raises(SlurmCommandNotFoundError, match="squeue"):
        SlurmClient().active_jobs(user="researcher")


@patch("OrthoEvol.Tools.slurm.client.shutil.which", return_value="/usr/bin/squeue")
@patch("OrthoEvol.Tools.slurm.client.subprocess.run")
def test_rejects_malformed_scheduler_output(
    mock_run: Mock,
    _mock_which: Mock,
) -> None:
    mock_run.return_value = subprocess.CompletedProcess(
        args=[],
        returncode=0,
        stdout="12345|orthologs|RUNNING\n",
        stderr="",
    )

    with pytest.raises(ValueError, match="Expected 7 Slurm fields"):
        SlurmClient().active_jobs(user="researcher")
