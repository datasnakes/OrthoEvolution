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


def archive(db_path, arch_path, config_file, delete=False):
    """
    Using YAML configuration, archive one or more directories recursively.

    This utility creates a YAML config dictionary that contains path-like objects for archiving.  The original data
    can be moved to the archive path or deleted all together.

    :param db_path:  A path to a folder that consists of the desired data.
    :param arch_path:  A path to an output folder for archived data.
    :param config_file:  The "Archive_Config" file.
    :param delete:  A flag for deleting the original data.  USE WITH CAUTION.
    :return:  Returns a list of paths to the *.tar.xz archive of the data and/or a path to the original data.
    """
    archive_list = []
    archive_log = LogIt().default(logname="Archive", logfile=None)
    with open(config_file, 'r') as yam:
        db_config_dict = yaml.safe_load(yam)

    archive_path = arch_path
    archive_dict = {}

    # Create a handle for creating separate archive files.
    if db_config_dict["Full"]:
        # For a full archive, individually archive each folder in the user database directory.
        full_path = db_path / archive_options["Full"]
        for folder in os.listdir(str(full_path)):
            archive_dict[folder] = db_path / Path(folder)
    else:
        # For custom archive options, set Full to False and select which data you would like to archive.
        for archive_key, archive_value in db_config_dict["Archive_Config"].items():
            # The desired projects for archiving must be explicitly listed in the config file.
            if archive_key == "Projects":
                if archive_value["flag"]:
                    for proj_key in archive_value.keys():
                        if proj_key != "flag":
                            if archive_value[proj_key]:
                                archive_dict[proj_key] = db_path / Path(proj_key)
            elif archive_value:
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
            archive_list.append(str(archive_filename))
        else:
            archive_log.critical("The original data will be moved recursively from %s to %s." % (data_path, output_pathname))
            output_pathname.mkdir()
            shutil.move(src=str(data_path), dst=str(output_pathname))
            shutil.move(src=str(archive_filename), dst=str(output_pathname))
            archive_list.append(str(output_pathname))

        Path(data_path).mkdir(parents=True, exist_ok=True)
    return archive_list


def get_size(start_path, units="KB"):
    """
    Determine the size of a directory or a file with the desired units.

    :param start_path:  A file or path for sizing up.
    :param units:  The denomination of bytes to return.
    :return:  The size as a string.  (e.g. "3.6 KB")
    """
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
