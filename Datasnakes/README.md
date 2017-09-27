Datasnakes Mini-Tutorial
==========================
Presently, this package is comprised of 4 major modules: Cookies, Manager, Orthologs,
and Tools. These 4 modules combine to offer a cohesive environment for easily creating,
managing, and deploying a data analysis pipeline of orthologous genes/species.

Additionally, this pipeline provides a way to easily visualize data, create/maintain project
structure, and edit structure/analysis tools to fit needs.

READMEs are provided in each module's directory, but we've compiled a mini tutorial here
that can inform users on how to use these modules.

Using Cookies
--------------
#### What is Cookies?
The Cookies module acts as a repository for custom [cookiecutter](https://github.com/audreyr/cookiecutter) templates.  

#### Why Cookies?
Each "Cookie" allows us to quickly template different parts of our project.  They are meant to help organize projects
and data in a standardized way.  The cookie module is used primarily by the Manager module.  

While these templates and the associated code in the Cookies/Manager modules are geared towards a CLI and a Flask framework, 
they are also highly useful for organizing multiple projects conducted by multiple users in a single repository.  This module is also useful for 
standalone projects.

#### Examples
**Templates used when creating a full repository:**
* _Cookies/new_repository_
* _Cookies/new_user_
* _Cookies/new_project_
* _Cookies/new_research_
* _Cookies/new_database_ (for NCBI, proprietary, etc. databases)
* _Cookies/new_app_ (for [R-Shiny](https://github.com/grabear/awesome-rshiny) applications)
* _Cookies/new_website_ (for [Flask](http://flask.pocoo.org/) applications)

**Template for standalone projects**
* _Cookies/new_basic_project_

Using Manager
--------------
#### What is Manager?

#### Why Manager?

#### Examples

Using Orthologs
----------------
#### What is Orthologs?

#### Why Orthologs?

#### Examples

Using Tools
------------
#### What is Tools?

#### Why Tools?

#### Examples