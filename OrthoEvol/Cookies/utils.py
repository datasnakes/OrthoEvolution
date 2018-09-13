import shutil
import datetime
import os
from pathlib import Path
from OrthoEvol.Tools.logit import LogIt


class CookieUtils(object):
    def __init__(self):
        self.archive_options = {
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

        self.bytesize_options = {
            "B": 1,
            "KB": 1024,
            "MB": 1048576,
            "GB": 1073741824,
            "TB": 1099511627776
        }

    def archive(self, database_path, archive_path, option, delete_flag=False):
        """
        Archive a database directory from a Cookie templated directory structure.

        This utility creates a YAML config dictionary that contains path-like
        objects for archiving.  The original data
        can be moved to the archive path or deleted all together.

        :param database_path:  A path to a folder that consists of the desired data.
        :param archive_path:  A path to an output folder for archived data.
        :param option:  An option for the archiving strategy.  Will be one of the keys in the archive_options.
        :param delete_flag:  A flag for deleting the original data.  USE WITH CAUTION.
        :return:  Returns a list of paths to the *.tar.xz archive of the data and/or a path to the original data.
        """
        archive_dict = {}
        archive_list = []
        archive_log = LogIt().default(logname="Archive", logfile=None)

        if option == "Full":
            full_path = Path(database_path) / self.archive_options["Full"]
            for folder in os.listdir(str(full_path)):
                if os.path.isdir(folder):
                    archive_dict[folder] = database_path / Path(folder)
        elif isinstance(option, list):
            for opt in option:
                other_path = Path(database_path) / self.archive_options[opt]
                archive_dict[opt] = other_path
        else:
            other_path = Path(database_path) / self.archive_options[option]
            archive_dict[option] = other_path

        for arch_name, data_path in archive_dict.items():
            root_dir = str(data_path.parent)
            base_dir = str(data_path.stem)
            d = datetime.datetime.now().strftime(format="%Y-%m-%d_%H%M")
            output_pathname = archive_path / Path(arch_name + "." + d)
            # Archive the desired data.
            data_size = self.get_size(start_path=str(data_path))
            archive_log.info("Archiving %s of data." % data_size)
            archive_filename = shutil.make_archive(base_name=str(output_pathname), format="xztar", root_dir=root_dir,
                                                   base_dir=base_dir, logger=archive_log)
            archive_size = self.get_size(archive_filename)
            archive_log.warning("A %s archive file was created at %s." % (archive_filename, archive_size))
            # TODO-ROB:  Logging.  And log to a README.md file.
            # Delete the files if desired.
            if delete_flag:
                archive_log.critical("The original data will be deleted recursively at %s." % data_path)
                from OrthoEvol import OrthoEvolWarning
                OrthoEvolWarning("You're about to delete your database (%s).  Are you sure??" % data_path)
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

    def get_size(self, start_path, units="KB"):
        """
        Determine the size of a directory or a file with the desired units.

        :param start_path:  A file or path for sizing up.
        :param units:  The denomination of bytes to return.
        :return:  The size as a string.  (e.g. "3.6 KB")
        """
        total_size = 0
        if os.path.isfile(start_path):
            size = os.path.getsize(start_path)
            size = str(size/self.bytesize_options[units]) + (" %s" % units)
            return size

        for dirpath, dirnames, filenames in os.walk(start_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total_size += os.path.getsize(fp)
        total_size = str(total_size/self.bytesize_options[units]) + (" %s" % units)
        return total_size
