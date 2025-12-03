"""Setup script for OrthoEvolution package.

This file is kept for programmatic setup tasks (e.g., BioSQL script permissions).
Most package metadata is now defined in pyproject.toml.

References:
https://github.com/pypa/sampleproject/blob/master/setup.py
https://github.com/biopython/biopython/blob/master/setup.py
http://python-packaging.readthedocs.io/en/latest/index.html
"""
from setuptools import setup
from setuptools.command.install import install
import os
from importlib import import_module
import pkg_resources


class PostInstallCommand(install):
    """Post-installation command to set BioSQL Perl script permissions."""

    def run(self):
        """Run the installation and then set script permissions."""
        install.run(self)
        # Set up the permissions for the BioSQL Perl scripts
        try:
            scripts = import_module("OrthoEvol.Manager.biosql.biosql_repo.scripts")
            biosql_scripts = pkg_resources.resource_filename(scripts.__name__, "")
            for file in os.listdir(biosql_scripts):
                if '.pl' in file:
                    script_path = os.path.join(biosql_scripts, file)
                    os.chmod(script_path, mode=755)
        except (ImportError, ModuleNotFoundError, FileNotFoundError):
            # Skip if module not available or files don't exist
            pass


# Most metadata is now in pyproject.toml
# This setup() call is minimal - setuptools will read pyproject.toml
setup(cmdclass={'install': PostInstallCommand})
