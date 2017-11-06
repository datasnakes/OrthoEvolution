import yaml
import shutil
import datetime
import os
from pathlib import Path
from Datasnakes.Tools import LogIt

archive_options = {
    "Full": Path(''),
    "NCBI": Path('NCBI'),
    "ITIS": Path('ITIS'),
    "NCBI_blast": Path('NCBI/blast'),
    "NCBI_blast_db": Path('NCBI/blast/db'),
    "NCBI_blast_windowmasker_files": Path('NCBI/blast/windowmasker_files'),
    "NCBI_pub_taxonomy": Path('NCBI/pub/taxonomy'),
    "NCBI_refseq_release": Path('NCBI/refseq/release'),
    "ITIS_taxonomy": Path('ITIS/taxonomy'),
}

bytesize_options = {
    "B": 1,
    "KB": 1024,
    "MB": 1048576,
    "GB": 1073741824,
    "TB": 1099511627776
}


def archive(db_path, config_file, delete=False):
    archive_log = LogIt().default(logname="Archive", logfile=None)
    with open(config_file, 'r') as yam:
        db_config_dict = yaml.safe_load(yam)

    archive_path = db_path / Path('archive')
    archive_dict = {}

    for archive_key, archive_value in db_config_dict["Archive_Config"].items():
        if archive_value:
            archive_dict[archive_key] = db_path / archive_options[archive_key]

    for arch_name, data_path in archive_dict.items():
        root_dir = str(data_path.parent)
        base_dir = str(data_path.stem)
        d = datetime.datetime.now().strftime(fmt="%Y-%m-%d_%H%M")
        output_pathname = archive_path / Path(arch_name + "." + d)
        # Archive the desired data.
        data_size = get_size(start_path=str(data_path))
        archive_log.info("Archiving %s of data." % data_size)
        archive_filename = shutil.make_archive(base_name=str(output_pathname), format="xztar", root_dir=root_dir,
                                               base_dir=base_dir)
        archive_size = get_size(archive_filename)
        archive_log.warning("A %s archive file was created at %s." % (archive_filename, archive_size))
        # TODO-ROB:  Logging.  And log to a README.md file.
        # Delete the files if desired.
        if delete:
            archive_log.critical("The original data will be deleted recursively at %s." % data_path)
            from Datasnakes import DatasnakesWarning
            DatasnakesWarning("You're about to delete your database (%s).  Are you sure??" % data_path)
            shutil.rmtree(path=data_path)
        else:
            archive_log.critical("The original data will be moved recursively from %s to %s." % (data_path, output_pathname))
            output_pathname.mkdir()
            shutil.move(src=str(data_path), dst=str(output_pathname))
            shutil.move(src=str(archive_filename), dst=str(output_pathname))

        Path(data_path).mkdir(parents=True, exist_ok=True)


def get_size(start_path, units="KB"):
    total_size = 0
    if os.path.isfile(start_path):
        size = os.path.getsize(start_path)
        size =str(size/bytesize_options[units]) + (" %s" % units)
        return size

    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    total_size = str(total_size/bytesize_options[units]) + (" %s" % units)
    return total_size
