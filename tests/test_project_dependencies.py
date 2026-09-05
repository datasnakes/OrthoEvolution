"""Tests for declared runtime dependencies."""

import tomllib
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_luigi_is_not_a_declared_dependency() -> None:
    pyproject_data = tomllib.loads(
        (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    dependencies = pyproject_data["project"]["dependencies"]
    requirements = (
        (REPOSITORY_ROOT / "requirements.txt")
        .read_text(encoding="utf-8")
        .splitlines()
    )

    assert all(not dependency.lower().startswith("luigi") for dependency in dependencies)
    assert all(not requirement.lower().startswith("luigi") for requirement in requirements)
