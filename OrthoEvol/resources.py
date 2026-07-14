"""Helpers for accessing files bundled with OrthoEvol packages."""

from importlib.resources import files
from pathlib import Path
from types import ModuleType


def package_resource_path(
    package: ModuleType,
    resource_name: str = "",
) -> Path:
    """Return a filesystem path for a resource in an installed package.

    OrthoEvol is configured as not zip-safe, so its installed resources are
    real files and can be passed to tools that require filesystem paths.
    """
    package_files = files(package)
    resource = (
        package_files.joinpath(resource_name) if resource_name else package_files
    )
    return Path(str(resource))
