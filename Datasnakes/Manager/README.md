Manager
=========
The classes and functions in this module have been designed to help manage
existing and new projects and pipelines using the Cookies module and the different utilities found in the Tools module.

### Management, RepoManagement, UserManagement, WebManagement, and Project Management classes

This module is intended to mesh with a Flask user interface.
* Whenever a new website is made the RepoManagement and WebManagement classes
are used.
    * Whenever a new user is created in the Flask webpage,
    the UserManagement class is used.
    * Whenever an existing user creates a new project,
    the ProjectManagement class is used.

However, this module does not have to be used to create a Flask
webpage.  The full repository can be used for higher level organization,
or standalone projects can be made using the ProjectManagements
_basic_project_ flag.

### DataManagement class

This module ties everything together into a pipeline.

#### **Cookies Use**

#### **Manager Use**

#### **Orthologs Use**

#### **Tools Use**

Examples
--------
In our [Examples module](https://github.com/datasnakes/Datasnakes-Scripts/tree/cookie_jar_patch/Examples),
you can see a perfect example of using Manager in **example_manager.py**.

```python
import os
from Datasnakes.Manager import DataManagement

DataManagement(pipeline="Ortho_CDS_1", start=True, new=True)
```