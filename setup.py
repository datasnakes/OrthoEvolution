"""Setup script for OrthoEvolution package.

This file is kept for programmatic setup tasks (e.g., BioSQL script permissions).
Most package metadata is now defined in pyproject.toml.

References:
https://github.com/pypa/sampleproject/blob/master/setup.py
https://github.com/biopython/biopython/blob/master/setup.py
http://python-packaging.readthedocs.io/en/latest/index.html
"""
from importlib import import_module
from pathlib import Path

from setuptools import setup
from setuptools.command.install import install


class PostInstallCommand(install):
    """Post-installation command to set BioSQL Perl script permissions."""

    def run(self) -> None:
        """Run the installation and then set script permissions."""
        install.run(self)
        # Set up the permissions for the BioSQL Perl scripts
        try:
            scripts = import_module("OrthoEvol.Manager.biosql.biosql_repo.scripts")
            biosql_scripts = Path(next(iter(scripts.__path__)))
            for script_path in biosql_scripts.glob("*.pl"):
                script_path.chmod(0o755)
        except (ImportError, ModuleNotFoundError, FileNotFoundError):
            # Skip if module not available or files don't exist
            pass


# Most metadata is now in pyproject.toml
# This setup() call is minimal - setuptools will read pyproject.toml
setup(cmdclass={'install': PostInstallCommand})
