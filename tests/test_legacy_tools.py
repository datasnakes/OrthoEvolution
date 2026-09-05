"""Minimal tests for small legacy command and integration helpers."""

from pathlib import Path
from unittest import mock
from zipfile import ZipFile

import pytest

from OrthoEvol.Manager.config.scripts.pipeline import OrthologPipeline
from OrthoEvol.Orthologs.command_line import main
from OrthoEvol.Tools.send2server.s2s import S2S
from OrthoEvol.Tools.slackify.notify import Slackify


def test_pipeline_iterates_over_gene_file(tmp_path: Path) -> None:
    index = tmp_path / "index"
    gene_path = tmp_path / "raw_data" / "GeneA"
    index.mkdir()
    gene_path.mkdir(parents=True)
    genes = tmp_path / "genes.csv"
    genes.write_text("GeneA\n", encoding="utf-8")
    (index / "pipeline.pbs").write_text("#!/bin/bash\n", encoding="utf-8")
    submitted: list[str] = []
    pipeline = OrthologPipeline(
        genes=genes,
        qsub_template="pipeline.pbs",
        worker_template="worker.py",
        home=index,
    )
    pipeline.submit = submitted.append

    pipeline.iterate()

    assert (gene_path / "GeneA.pbs").is_file()
    assert len(submitted) == 1
    assert "GENE=GeneA" in submitted[0]
    assert str(index / "worker.py") in submitted[0]


def test_placeholder_command_line(capsys: pytest.CaptureFixture[str]) -> None:
    main()

    assert capsys.readouterr().out == "Placeholder\n"


def test_s2s_builds_commands_and_zips_local_files(tmp_path: Path) -> None:
    source = tmp_path / "input.txt"
    source.write_text("sequence data", encoding="utf-8")
    sender = S2S(
        username="user",
        server_address="example.org",
        dest_path="/data",
        comp_filename="archive.zip",
        zip_path=tmp_path,
        compressed=False,
    )
    sender.ignore_parts = []

    archive = Path(sender.to_zip())

    assert sender.send_cmd == "scp archive.zip user@example.org:/data"
    with ZipFile(archive) as zip_file:
        assert zip_file.read("input.txt").decode() == "sequence data"


@pytest.mark.parametrize(("method_name", "status", "message"), [
    ("scpto", 0, "file sent"),
    ("cpto", 1, "file not sent"),
])
def test_s2s_reports_transfer_status(
    method_name: str,
    status: int,
    message: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    sender = S2S(
        dest_path=tmp_path,
        comp_filename="archive.zip",
        zip_path=tmp_path,
        compressed=False,
    )
    monkeypatch.setattr(
        "OrthoEvol.Tools.send2server.s2s.subprocess.call",
        lambda *args, **kwargs: status,
    )

    getattr(sender, method_name)("archive.zip")

    assert message in capsys.readouterr().out


def test_slackify_reads_config_and_delegates_helpers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = tmp_path / "slack.cfg"
    config.write_text("[APIKEYS]\nslack = token\n", encoding="utf-8")
    slack = mock.Mock()
    slack.channels.get_channel_id.return_value = "C123"
    slack.users.list.return_value.body = {"members": [{"name": "alice"}]}
    slack.channels.list.return_value.body = {"channels": [{"name": "general"}]}
    monkeypatch.setattr(
        "OrthoEvol.Tools.slackify.notify.Slacker", lambda api_key: slack
    )
    notifier = Slackify(config)

    notifier.upload_file("report.txt", "general")
    notifier.send_msg("general", "finished")

    slack.files.upload.assert_called_once_with(
        file_="report.txt", channels="C123"
    )
    slack.chat.post_message.assert_called_once_with(
        "general", "finished", as_user=True
    )
    assert notifier.list_users() == ["alice"]
    assert notifier.list_channels() == ["general"]


def test_slackify_rejects_missing_config(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="configuriation file not found"):
        Slackify(tmp_path / "missing.cfg")
