#!/usr/bin/env python3
# ---
# [Shaurita Hutchins] update_docs.py
# Automated documentation generation and update script for OrthoEvolution.
# ---

"""Automated documentation generation and update script.

This script automates the process of:
1. Reading version from setup.py
2. Updating Sphinx configuration
3. Converting README.md files to RST (optional)
4. Regenerating API documentation using sphinx-apidoc
5. Building documentation
"""

import os
import re
import subprocess
import sys
import configparser
from pathlib import Path
from argparse import ArgumentParser
from typing import List, Dict


def get_version_from_setup():
    """Extract version from setup.py.

    :return: Version string from setup.py
    """
    setup_path = Path(__file__).parent.parent / 'setup.py'
    with open(setup_path, 'r') as f:
        content = f.read()
        match = re.search(r"version\s*=\s*['\"]([^'\"]+)['\"]", content)
        if match:
            return match.group(1)
    raise ValueError("Could not find version in setup.py")


def update_conf_py():
    """Update conf.py with current version (already automated in conf.py)."""
    print("✓ conf.py automatically reads version from setup.py")


def convert_readmes():
    """Convert README.md files to RST format."""
    try:
        import pypandoc
    except ImportError:
        print("\n⚠ Warning: pypandoc not installed.")
        print("  Skipping README.md conversion. Install with: pip install pypandoc")
        return

    docs_dir = Path(__file__).parent
    package_dir = docs_dir.parent / 'OrthoEvol'
    docs_source = docs_dir / 'docs' / 'source'
    config_file = docs_dir / 'docs.cfg'

    print("\n" + "-" * 60)
    print("Converting README.md files to RST...")

    # Skip patterns for cookiecutter templates (exact directory name matches)
    skip_dirs = {
        'new_basic_project', 'new_repository', 'new_database_repo',
        'new_research', 'new_user', 'new_website',
        '__pycache__', '.git'
    }
    
    # Skip if path contains these patterns
    skip_path_patterns = {
        '{{cookiecutter',
        '__init__.py'
    }

    # Read docs.cfg mapping
    file_mapping = {}
    if config_file.exists():
        config = configparser.ConfigParser()
        config.read(config_file)
        for section in config.sections():
            file_mapping[section] = [
                config[section][key] for key in config[section]
            ]

    # Create reverse mapping: filename -> target directory
    filename_to_dir = {}
    for section, filenames in file_mapping.items():
        for filename in filenames:
            filename_to_dir[filename] = section

    # Find README.md files
    if not package_dir.exists():
        print(f"  Error: Package directory {package_dir} not found.")
        return

    files2convert = []
    for readme_file in package_dir.rglob('README.md'):
        # Skip if any parent directory name is in skip list
        path_parts = readme_file.parts
        if any(part in skip_dirs for part in path_parts):
            continue
        # Skip if path contains skip patterns
        readme_str = str(readme_file)
        if any(pattern in readme_str for pattern in skip_path_patterns):
            continue
        files2convert.append(readme_file)

    if not files2convert:
        print("  No README.md files found to convert.")
        return

    print(f"  Found {len(files2convert)} README.md files to process...")

    converted_count = 0
    skipped_count = 0

    for file2convert in sorted(files2convert):
        parent_dir = file2convert.parent.name.lower()
        expected_filename = f"{parent_dir}readme.rst"

        # Determine output directory
        if expected_filename in filename_to_dir:
            section = filename_to_dir[expected_filename]
            outpath = docs_source / section
        elif (parent_dir == 'orthologs' and
              file2convert.parent.parent.name.lower() == 'orthoevol'):
            outpath = docs_source / 'orthologs'
        elif (parent_dir == 'tools' and
              file2convert.parent.parent.name.lower() == 'orthoevol'):
            outpath = docs_source / 'tools'
        elif (file2convert.name == 'README.md' and
              file2convert.parent.name == 'OrthoEvol'):
            outpath = docs_source / 'tutorial'
            expected_filename = 'orthoevolreadme.rst'
        else:
            skipped_count += 1
            continue

        # Convert file
        outpath.mkdir(parents=True, exist_ok=True)
        final_filename = expected_filename
        outfile = outpath / final_filename

        try:
            pypandoc.convert_file(
                str(file2convert),
                'rst',
                outputfile=str(outfile)
            )
            print(f"  ✓ Converted: {file2convert.name} -> {outfile.name}")
            converted_count += 1
        except Exception as e:
            print(f"  ✗ Error converting {file2convert}: {e}")

    print(f"✓ Conversion complete: {converted_count} files converted, {skipped_count} skipped.")


def regenerate_api_docs():
    """Regenerate API documentation using sphinx-apidoc.

    This removes old module docs and regenerates them from the source code.
    """
    docs_dir = Path(__file__).parent
    source_dir = docs_dir / 'docs' / 'source'
    modules_dir = source_dir / 'modules'
    package_dir = docs_dir.parent / 'OrthoEvol'

    print(f"Regenerating API documentation...")
    print(f"  Source: {package_dir}")
    print(f"  Output: {modules_dir}")

    # Remove old module files
    if modules_dir.exists():
        for old_file in modules_dir.glob('*.rst'):
            if old_file.name != 'modules.rst':
                old_file.unlink()
                print(f"  Removed: {old_file.name}")

    # Run sphinx-apidoc
    cmd = [
        'sphinx-apidoc',
        '-o', str(modules_dir),
        str(package_dir),
        '--separate',
        '--force'
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✓ API documentation regenerated successfully")
        if result.stdout:
            print(f"  Output: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error regenerating API docs: {e}")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("✗ sphinx-apidoc not found. Install sphinx: pip install sphinx")
        sys.exit(1)


def build_docs(builder='html'):
    """Build the documentation.

    :param builder: Sphinx builder to use (default: html)
    """
    docs_dir = Path(__file__).parent
    source_dir = docs_dir / 'docs' / 'source'
    build_dir = docs_dir / 'docs' / '_build'

    print(f"\nBuilding documentation ({builder})...")
    print(f"  Source: {source_dir}")
    print(f"  Build: {build_dir}")

    cmd = ['sphinx-build', '-b', builder, str(source_dir), str(build_dir)]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"✓ Documentation built successfully")
        print(f"  Output directory: {build_dir}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error building docs: {e}")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("✗ sphinx-build not found. Install sphinx: pip install sphinx")
        sys.exit(1)


def main():
    """Main function to run all documentation updates."""
    parser = ArgumentParser(description='Update OrthoEvolution documentation')
    parser.add_argument(
        '--skip-readmes',
        action='store_true',
        help='Skip README.md to RST conversion'
    )
    parser.add_argument(
        '--skip-build',
        action='store_true',
        help='Skip building documentation (only regenerate)'
    )
    args = parser.parse_args()

    print("=" * 60)
    print("OrthoEvolution Documentation Update Script")
    print("=" * 60)

    # Get current version
    try:
        version = get_version_from_setup()
        print(f"\nCurrent version: {version}")
    except ValueError as e:
        print(f"✗ {e}")
        sys.exit(1)

    # Update configuration (already automated)
    update_conf_py()

    # Convert README.md files (optional)
    if not args.skip_readmes:
        convert_readmes()

    # Regenerate API docs
    print("\n" + "-" * 60)
    regenerate_api_docs()

    # Build docs
    if not args.skip_build:
        print("\n" + "-" * 60)
        build_docs()

    print("\n" + "=" * 60)
    print("Documentation update complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()
