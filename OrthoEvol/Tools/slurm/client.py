"""Small, synchronous wrappers around standard Slurm commands."""

import getpass
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Final

SLURM_FIELD_DELIMITER: Final = "|"
SLURM_OUTPUT_FIELDS: Final = (
    "job_id",
    "name",
    "state",
    "partition",
    "elapsed",
    "nodes",
    "node_list_or_reason",
)


class SlurmCommandNotFoundError(RuntimeError):
    """Report that a required Slurm executable is unavailable."""


@dataclass(frozen=True, slots=True)
class SlurmJob:
    """Normalized job information shared by active and historical queries."""

    job_id: str
    name: str
    state: str
    partition: str
    elapsed: str
    nodes: int | None
    node_list_or_reason: str


def _parse_node_count(raw_node_count: str) -> int | None:
    """Keep unavailable node counts distinct from a real allocation of zero."""
    if not raw_node_count:
        return None

    try:
        return int(raw_node_count)
    except ValueError as error:
        raise ValueError(
            f"Invalid Slurm node count: {raw_node_count!r}"
        ) from error


def _parse_slurm_rows(output: str) -> list[SlurmJob]:
    """Parse the seven-field contract requested from squeue or sacct."""
    jobs: list[SlurmJob] = []
    expected_field_count = len(SLURM_OUTPUT_FIELDS)

    for line_number, line in enumerate(output.splitlines(), start=1):
        if not line.strip():
            continue

        fields = [
            value.strip() for value in line.split(SLURM_FIELD_DELIMITER)
        ]
        if len(fields) != expected_field_count:
            raise ValueError(
                f"Expected {expected_field_count} Slurm fields on line "
                f"{line_number}, found {len(fields)}: {line!r}"
            )

        jobs.append(
            SlurmJob(
                job_id=fields[0],
                name=fields[1],
                state=fields[2],
                partition=fields[3],
                elapsed=fields[4],
                nodes=_parse_node_count(fields[5]),
                node_list_or_reason=fields[6],
            )
        )

    return jobs


class SlurmClient:
    """Run one-shot Slurm submission and inspection commands."""

    _squeue_format: Final = "%i|%j|%T|%P|%M|%D|%R"
    _sacct_format: Final = (
        "JobIDRaw,JobName,State,Partition,Elapsed,AllocNodes,NodeList"
    )

    @staticmethod
    def _run(command: list[str]) -> str:
        """Execute a Slurm command without invoking a shell."""
        executable = command[0]
        if shutil.which(executable) is None:
            raise SlurmCommandNotFoundError(
                f"Required Slurm command {executable!r} was not found in PATH."
            )

        completed_process = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed_process.stdout

    def submit(self, script: Path) -> str:
        """Submit an existing batch script and return its Slurm job ID."""
        script_path = Path(script)
        if not script_path.is_file():
            raise FileNotFoundError(
                f"Slurm batch script does not exist: {script_path}"
            )

        output = self._run(["sbatch", "--parsable", str(script_path)])
        job_id = output.strip().partition(";")[0]
        if not job_id:
            raise ValueError("sbatch returned an empty job ID.")
        return job_id

    def active_jobs(self, user: str | None = None) -> list[SlurmJob]:
        """Return one snapshot of active jobs for a single user."""
        username = user or getpass.getuser()
        output = self._run(
            [
                "squeue",
                "--noheader",
                f"--user={username}",
                f"--format={self._squeue_format}",
            ]
        )
        return _parse_slurm_rows(output)

    def job_history(self, job_id: str) -> list[SlurmJob]:
        """Return the allocation-level accounting record for one job."""
        if not job_id.strip():
            raise ValueError("A Slurm job ID is required.")

        output = self._run(
            [
                "sacct",
                "--allocations",
                "--noheader",
                "--parsable2",
                f"--jobs={job_id}",
                f"--format={self._sacct_format}",
            ]
        )
        return _parse_slurm_rows(output)
