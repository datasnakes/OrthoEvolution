"""Tests for classes that consume bundled OrthoEvol resources."""

import os
import tempfile
import unittest
from pathlib import Path
from shutil import rmtree
from types import SimpleNamespace
from unittest import mock

from OrthoEvol.Manager.config import yml
from OrthoEvol.Manager.biosql.biosql import BaseBioSQL
from OrthoEvol.Manager.biosql.biosql_repo import scripts as sql_scripts
from OrthoEvol.Manager.data_management import DataMana
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Orthologs.Blast.comparative_genetics import (
    BaseComparativeGenetics,
)
from OrthoEvol.Orthologs.Phylogenetics.PAML.codeml import CodemlRun
from OrthoEvol.Tools.pbs.qstat import BaseQstat
from OrthoEvol.Tools.pbs.qsub import Qsub


class TestResourceConsumers(unittest.TestCase):
    """Verify callers resolve the bundled files they depend on."""

    def setUp(self) -> None:
        """Create an isolated filesystem for constructor smoke tests."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.addCleanup(rmtree, self.test_dir, ignore_errors=True)

    @mock.patch("OrthoEvol.Manager.data_management.package_resource_path")
    def test_data_manager_selects_pipeline_config(
        self,
        mock_resource_path: mock.Mock,
    ) -> None:
        """Select the correct bundled configuration for new and existing runs."""
        expected_calls = (
            (True, "pipeline.yml"),
            (False, "initialize_old.yml"),
        )

        for is_new, expected_name in expected_calls:
            with self.subTest(is_new=is_new):
                DataMana(pipeline="Ortho_CDS_1", new=is_new)

                mock_resource_path.assert_called_with(yml, expected_name)

    @mock.patch("OrthoEvol.Manager.biosql.biosql.package_resource_path")
    def test_biosql_resolves_taxonomy_scripts(
        self,
        mock_resource_path: mock.Mock,
    ) -> None:
        """Resolve the BioSQL script directory and both taxonomy loaders."""
        mock_resource_path.side_effect = (
            self.test_dir,
            self.test_dir / "load_ncbi_taxonomy.pl",
            self.test_dir / "load_itis_taxonomy.pl",
        )

        biosql = BaseBioSQL(
            database_name="test.db",
            project="test-project",
            project_path=self.test_dir,
            proj_mana=None,
        )

        self.assertEqual(biosql.scripts, self.test_dir)
        self.assertEqual(
            biosql.ncbi_taxon_script,
            self.test_dir / "load_ncbi_taxonomy.pl",
        )
        self.assertEqual(
            biosql.itis_taxon_script,
            self.test_dir / "load_itis_taxonomy.pl",
        )
        self.assertEqual(
            mock_resource_path.call_args_list,
            [
                mock.call(sql_scripts),
                mock.call(sql_scripts, "load_ncbi_taxonomy.pl"),
                mock.call(sql_scripts, "load_itis_taxonomy.pl"),
            ],
        )

    def test_scheduler_classes_resolve_configuration_files(self) -> None:
        """Resolve scheduler configuration without invoking cluster commands."""
        qstat_config = self.test_dir / "qstat.yml"
        pbs_template = self.test_dir / "temp.pbs"

        with mock.patch(
            "OrthoEvol.Tools.pbs.qstat.package_resource_path",
            return_value=qstat_config,
        ) as mock_qstat_resource:
            qstat = BaseQstat(job_id="123.server", home=self.test_dir / "qstat")

        with mock.patch(
            "OrthoEvol.Tools.pbs.qsub.package_resource_path",
            return_value=pbs_template,
        ) as mock_qsub_resource:
            qsub = Qsub(
                job_name="test-job",
                base_job_id="abcde",
                pbs_working_dir=self.test_dir,
                python_script="analysis.py",
                pbs_command_list=[],
            )

        self.assertEqual(qstat._yaml_config, qstat_config)
        self.assertEqual(qsub.pbs_template, pbs_template)
        mock_qstat_resource.assert_called_once()
        mock_qsub_resource.assert_called_once()

    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.shutil.copy"
    )
    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.package_resource_path"
    )
    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.FullUtilities"
    )
    def test_comparative_genetics_resolves_packaged_accessions(
        self,
        mock_utilities: mock.Mock,
        mock_resource_path: mock.Mock,
        mock_copy: mock.Mock,
    ) -> None:
        """Skip packaged accession copying when no accession file is configured."""
        accession_path = self.test_dir / "accessions.csv"
        mock_resource_path.return_value = accession_path
        mock_utilities.return_value.attribute_config.return_value = SimpleNamespace(
            project_index=self.test_dir
        )

        analysis = BaseComparativeGenetics(
            project="test-project",
            project_path=self.test_dir,
            acc_file=None,
            proj_mana=None,
            copy_from_package=True,
        )

        mock_resource_path.assert_not_called()
        mock_copy.assert_not_called()
        self.assertEqual(analysis.project, "test-project")
        self.assertEqual(analysis.project_path, self.test_dir / "test-project")

    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.shutil.copy"
    )
    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.package_resource_path"
    )
    @mock.patch(
        "OrthoEvol.Orthologs.Blast.comparative_genetics.FullUtilities"
    )
    def test_comparative_genetics_derives_project_from_existing_path(
        self,
        mock_utilities: mock.Mock,
        mock_resource_path: mock.Mock,
        mock_copy: mock.Mock,
    ) -> None:
        """Use the final path component as the project name."""
        existing_project_path = self.test_dir / "existing-project"
        accession_path = self.test_dir / "accessions.csv"
        mock_resource_path.return_value = accession_path
        mock_utilities.return_value.attribute_config.return_value = SimpleNamespace(
            project_index=self.test_dir
        )

        analysis = BaseComparativeGenetics(
            project=None,
            project_path=existing_project_path,
            acc_file=None,
            proj_mana=None,
            copy_from_package=True,
        )

        self.assertEqual(analysis.project, "existing-project")
        self.assertEqual(analysis.project_path, existing_project_path)
        mock_utilities.return_value.attribute_config.assert_called_once_with(
            cls=analysis,
            composer=None,
            checker=ProjectManagement,
            project="existing-project",
            project_path=self.test_dir,
        )
        mock_resource_path.assert_not_called()
        mock_copy.assert_not_called()

    @mock.patch(
        "OrthoEvol.Orthologs.Phylogenetics.PAML.codeml.package_resource_path"
    )
    @mock.patch("OrthoEvol.Orthologs.Phylogenetics.PAML.codeml.copy")
    @mock.patch("OrthoEvol.Orthologs.Phylogenetics.PAML.codeml.codeml.Codeml")
    def test_codeml_resolves_control_template(
        self,
        mock_codeml: mock.Mock,
        mock_copy: mock.Mock,
        mock_resource_path: mock.Mock,
    ) -> None:
        """Resolve the bundled control template before configuring Codeml."""
        original_directory = Path.cwd()
        self.addCleanup(os.chdir, original_directory)
        control_template = self.test_dir / "codeml.ctl"
        mock_resource_path.return_value = control_template
        mock_copy.side_effect = lambda source, destination: str(
            Path(destination) / Path(source).name
        )

        codeml_run = CodemlRun(
            P2N_alignment="alignment.phy",
            iqtree_newick="gene_iqtree.nwk",
            home=self.test_dir,
        )

        self.assertEqual(codeml_run.control_template, control_template)
        mock_resource_path.assert_called_once()
        mock_codeml.return_value.read_ctl_file.assert_called_once_with(
            control_template
        )


if __name__ == "__main__":
    unittest.main()
