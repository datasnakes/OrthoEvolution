from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher
from OrthoEvol.Manager.config import yml
from pkg_resources import resource_filename
import yaml

db_config_file = resource_filename(yml.__name__, "database_config.yml")
pm_config_file = resource_filename(yml.__name__, "config_template_existing.yml")
with open(pm_config_file, 'r') as f:
    pm_config = yaml.safe_load(f)
pm = ProjectManagement(**pm_config["Management_config"])

dd = DatabaseDispatcher(db_config_file, pm)

dd.dispatch(dd.strategies, dd.dispatcher, dd.configuration)
dd.uploading_dispatch(["NCBI_refseq_release"])
