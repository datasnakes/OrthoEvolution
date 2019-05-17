# Manager Documentation

The classes and functions in this module have been designed to help manage
existing and new projects using the Cookies module as well as
the different utilities found in the Tools module.

## Why a manager?

This module is intended to mesh with a Flask user interface.  While the 
Flask server/client interface is not currently set up, we are developing
the package so that it will be easier to implement.

* Whenever a new website is made the RepoManagement and WebManagement classes
are used.
  * Whenever a new user is created in the Flask webpage,
    the UserManagement class is used.
  * Whenever an existing user creates a new project,
    the ProjectManagement class is used.

This module does not have to be used to create a Flask
webpage.  The full repository can be used for higher level organization,
or standalone projects can be made using the ProjectManagements
_basic_project_ flag.

## Example

Using DatabaseDispatcher to set up the proper databases using a `YAML`
config file:

```python
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher
from OrthoEvol.Manager.config import yml
from pkg_resources import resource_filename
from pathlib import Path
import yaml

# Set up project management
pm_config_file = resource_filename(yml.__name__, "initialize_new.yml")
with open(pm_config_file, 'r') as f:
    pm_config = yaml.load(f)
pm = ProjectManagement(**pm_config["Management_config"])

# Set up database management
db_config_file = resource_filename(yml.__name__, "databases.yml")
with open(db_config_file, 'r') as f:
    db_config = yaml.load(f)
config = db_config.update(pm_config)

# Generate main config file for this job
config_file = pm.user_log / Path("upload_config.yml")
with open(str(config_file), 'w') as cf:
    yaml.dump(config, cf, default_flow_style=False)

# Set up database dispatcher and dispatch the functions
dd = DatabaseDispatcher(config_file, pm)
dd.dispatch(dd.strategies, dd.dispatcher, dd.configuration)
```

## Notes

Please view our [BioSQL documentation](https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager/BioSQL/README.md) and view some of the
static/config related [files](https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager/config/).