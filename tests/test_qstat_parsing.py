"""Regression tests for parsing OpenPBS qstat text output."""

import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path
from shutil import rmtree

from OrthoEvol.Tools.pbs.qstat import BaseQstat


class TestQstatKeywordParsing(unittest.TestCase):
    """Verify qstat fields and wrapped values retain their legacy shape."""

    def setUp(self) -> None:
        """Create a parser whose output directory is isolated per test."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.addCleanup(rmtree, self.test_dir, ignore_errors=True)
        self.qstat = BaseQstat(
            job_id="4.scheduler",
            home=self.test_dir,
        )

    def test_parses_fields_and_wrapped_variable_list(self) -> None:
        """Distinguish four-space fields from deeper continuation lines."""
        variable_start = "    Variable_List = PBS_O_HOME=/home/adminuser,\n"
        variable_continuation = "        PBS_O_LOGNAME=adminuser,PBS_O_QUEUE=workq\n"
        job_data = OrderedDict(
            {
                "4.scheduler": [
                    "    Job_Name = STDIN\n",
                    "    Job_Owner = adminuser@submit.example.org\n",
                    "    resources_used.cpupercent = 0\n",
                    "    job_state = R\n",
                    "    Resource_List.ncpus = 1\n",
                    variable_start,
                    variable_continuation,
                    "\n",
                ]
            }
        )

        parsed_jobs = self.qstat.identify_qstat_keywords(job_data)

        self.assertEqual(parsed_jobs["4.scheduler"]["job_state"], "    job_state = R\n")
        self.assertEqual(
            parsed_jobs["4.scheduler"]["Variable_List"],
            [variable_start, variable_continuation],
        )

    def test_rejects_unknown_attribute(self) -> None:
        """Report new scheduler fields instead of silently dropping them."""
        job_data = {"4.scheduler": ["    custom_attribute = value\n"]}

        with self.assertRaisesRegex(KeyError, "custom_attribute"):
            self.qstat.identify_qstat_keywords(job_data)

    def test_accepts_explicit_extra_attribute(self) -> None:
        """Allow callers to opt into site-specific qstat attributes."""
        raw_line = "    custom_attribute = value\n"
        job_data = {"4.scheduler": [raw_line]}

        parsed_jobs = self.qstat.identify_qstat_keywords(
            job_data,
            extra_keywords=["custom_attribute"],
        )

        self.assertEqual(
            parsed_jobs["4.scheduler"]["custom_attribute"],
            raw_line,
        )

    def test_rejects_continuation_without_an_attribute(self) -> None:
        """Reject wrapped content that has no preceding field to own it."""
        job_data = {"4.scheduler": ["        orphaned content\n"]}

        with self.assertRaisesRegex(ValueError, "orphaned content"):
            self.qstat.identify_qstat_keywords(job_data)

    def test_rejects_non_mapping_input(self) -> None:
        """Fail clearly when Stage 1 output is not supplied."""
        with self.assertRaisesRegex(TypeError, "dictionary"):
            self.qstat.identify_qstat_keywords([])  # type: ignore[arg-type]

    def test_complete_parser_converts_nested_values(self) -> None:
        """Exercise the parsing stages with one representative scheduler job."""
        raw_data = [
            "Job Id: 4.scheduler\n",
            "    Job_Name = analysis\n",
            "    Resource_List.ncpus = 4\n",
            "    Variable_List = PBS_O_HOME=/home/user,\n",
            "        PBS_O_QUEUE=workq\n",
        ]

        parsed = self.qstat.to_dict(raw_data, ordered=False)

        self.assertEqual(parsed["4.scheduler"]["Job_Name"], "analysis")
        self.assertEqual(parsed["4.scheduler"]["Resource_List"]["ncpus"], 4)
        self.assertEqual(
            parsed["4.scheduler"]["Variable_List"]["PBS_O_QUEUE"], "workq"
        )


if __name__ == "__main__":
    unittest.main()
