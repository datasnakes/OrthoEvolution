"""Regression tests for parsing OpenPBS qstat text output."""

import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path
from shutil import rmtree
from types import SimpleNamespace
from unittest import mock

from OrthoEvol.Tools.pbs.qstat import BaseQstat, MultiQstat, Qstat


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

    def test_target_and_dataframe_filter_requested_job(self) -> None:
        """Keep target selection and dynamic dataframe generation aligned."""
        jobs = {
            "4.scheduler": {
                "Job_Name": "analysis",
                "job_state": "R",
                "resources_used.cpupercent": 75,
            }
        }

        target = self.qstat.target_data(jobs, "4.scheduler")
        dataframe = self.qstat.to_dataframe(jobs, "4.scheduler")

        self.assertEqual(target["Job_Name"], "analysis")
        self.assertEqual(dataframe.loc[0, "Job_Id"], "4.scheduler")
        self.assertEqual(dataframe.loc[0, "resources_used.cpupercent"], 75)

    def test_target_data_rejects_missing_job(self) -> None:
        """Report a missing target before downstream formatting."""
        with self.assertRaisesRegex(ValueError, "target job does not exist"):
            self.qstat.target_data({}, "missing.scheduler")

    def test_static_data_parses_scheduler_timestamps(self) -> None:
        """Normalize timestamps while excluding changing resource usage."""
        jobs = {
            "4.scheduler": {
                "Job_Name": "analysis",
                "ctime": "Thu Sep 4 12:00:00 2026",
                "resources_used.cpupercent": 75,
            }
        }

        static = self.qstat.static_data(jobs, "4.scheduler")

        self.assertEqual(static["4.scheduler"]["Job_Name"], "analysis")
        self.assertIn("2026-09-04 12:00:00", static["4.scheduler"]["ctime"])
        self.assertNotIn("resources_used.cpupercent", static["4.scheduler"])

    def test_configure_data_file_preserves_one_header(self) -> None:
        """Merge CSV history without repeating its header."""
        source = self.test_dir / "source.csv"
        output = self.test_dir / "combined.csv"
        source.write_text("job,state\n4,R\n", encoding="utf-8")

        self.qstat.configure_data_file(output, source)
        self.qstat.configure_data_file(output, source)

        self.assertEqual(
            output.read_text(encoding="utf-8").splitlines(),
            ["job,state", "4,R", "4,R"],
        )

    def test_qstat_output_reads_text_and_json(self) -> None:
        """Read both output formats after the command boundary writes them."""
        def system_cmd(*args: object, **kwargs: object) -> SimpleNamespace:
            log_file = Path(str(kwargs["file_name"]))
            content = str(kwargs.pop("content")) if "content" in kwargs else ""
            log_file.write_text(content, encoding="utf-8")
            return SimpleNamespace(returncode=0)

        for capture_json, content, expected in (
            (False, "Job Id: 4.scheduler\n", ["Job Id: 4.scheduler\n"]),
            (True, '{"Jobs": {"4.scheduler": {}}}', {"Jobs": {"4.scheduler": {}}}),
        ):
            self.qstat.qstat_utils.system_cmd = (
                lambda *args, _content=content, **kwargs: system_cmd(
                    *args, content=_content, **kwargs
                )
            )
            result = self.qstat.qstat_output(
                self.qstat.cmd,
                self.qstat.qstat_log_file,
                capture_json=capture_json,
            )
            self.assertEqual(result, expected)

    def test_run_qstat_writes_local_reports(self) -> None:
        """Coordinate JSON parsing and local report generation."""
        jobs = {
            "Jobs": {
                "4.scheduler": {
                    "Job_Name": "analysis",
                    "job_state": "R",
                }
            }
        }
        self.qstat.qstat_output = lambda **kwargs: jobs

        self.qstat.run_qstat(capture_json=True)

        self.assertEqual(self.qstat.job_dict["Job_Name"], "analysis")
        self.assertTrue(self.qstat.data_file.is_file())
        self.assertTrue(self.qstat.info_file.is_file())

    def test_local_reports_can_be_overwritten(self) -> None:
        """Replace existing CSV and YAML reports when explicitly requested."""
        jobs = {"4.scheduler": {"Job_Name": "analysis", "job_state": "R"}}
        csv_file = self.test_dir / "job.csv"
        yaml_file = self.test_dir / "job.yml"
        csv_file.write_text("old data\n", encoding="utf-8")
        yaml_file.write_text("old data\n", encoding="utf-8")

        self.qstat.to_csv(csv_file, jobs, "4.scheduler", overwrite=True)
        self.qstat.static_data_to_yaml(
            yaml_file, jobs, "4.scheduler", overwrite=True
        )

        self.assertNotIn("old data", csv_file.read_text(encoding="utf-8"))
        self.assertNotIn("old data", yaml_file.read_text(encoding="utf-8"))

    def test_constructor_supports_json_and_named_data_files(self) -> None:
        """Configure common command and output-path variants."""
        json_qstat = BaseQstat(
            job_id="4.scheduler",
            home=self.test_dir / "json",
            capture_json=True,
        )
        input_qstat = BaseQstat(
            job_id="4.scheduler",
            home=self.test_dir / "input",
            infile="history.csv",
        )
        output_qstat = BaseQstat(
            job_id="4.scheduler",
            home=self.test_dir / "output",
            outfile="results.csv",
        )

        self.assertEqual(json_qstat.cmd, "qstat -f 4 -F json")
        self.assertEqual(input_qstat.data_file.name, "history.csv")
        self.assertEqual(output_qstat.data_file.name, "results.csv")

    def test_local_reports_append_to_existing_files(self) -> None:
        """Append new observations when overwrite is disabled."""
        jobs = {"4.scheduler": {"Job_Name": "analysis", "job_state": "R"}}
        csv_file = self.test_dir / "append.csv"
        yaml_file = self.test_dir / "append.yml"
        csv_file.write_text("existing\n", encoding="utf-8")
        yaml_file.write_text("existing\n", encoding="utf-8")

        self.qstat.to_csv(csv_file, jobs, "4.scheduler")
        self.qstat.static_data_to_yaml(yaml_file, jobs, "4.scheduler")

        self.assertTrue(csv_file.read_text(encoding="utf-8").startswith("existing\n"))
        self.assertTrue(yaml_file.read_text(encoding="utf-8").startswith("existing\n"))

    def test_run_qstat_supports_text_and_sqlite_dispatch(self) -> None:
        """Parse text output and invoke the optional SQLite stage."""
        self.qstat.qstat_output = lambda **kwargs: [
            "Job Id: 4.scheduler\n",
            "    Job_Name = analysis\n",
            "    job_state = R\n",
        ]
        self.qstat.to_sqlite = mock.Mock()

        self.qstat.run_qstat(csv_flag=False, sqlite_flag=True)

        self.assertEqual(self.qstat.job_dict["Job_Name"], "analysis")
        self.qstat.to_sqlite.assert_called_once_with()

    def test_qstat_countdown_and_watch_delegate_without_waiting(self) -> None:
        """Exercise short control methods without polling a scheduler."""
        qstat = Qstat(job_id="4.scheduler", home=self.test_dir, wait_time=1)
        qstat._watch = mock.Mock()

        with mock.patch("OrthoEvol.Tools.pbs.qstat.sleep") as sleep:
            qstat.countdown(1)
        qstat.watch(max_count=1)

        sleep.assert_called_once_with(1)
        qstat._watch.assert_called_once_with(count=None, max_count=1)

    def test_multi_qstat_builds_one_watcher_per_job(self) -> None:
        """Create watcher objects without starting asynchronous polling."""
        multi_qstat = MultiQstat(
            jobs=["4.scheduler", "5.scheduler"], config_home=self.test_dir
        )

        watchers = multi_qstat.get_qstat_dict(
            ["4.scheduler", "5.scheduler"], wait_time=10
        )

        self.assertEqual(set(watchers), {"4.scheduler", "5.scheduler"})
        self.assertTrue(all(watcher.wait_time == 10 for watcher in watchers.values()))


if __name__ == "__main__":
    unittest.main()
