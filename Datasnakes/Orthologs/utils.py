from pathlib import Path
import os
from Datasnakes.Cookies import Oven

def config_composition(instance, project_path, project, composer, checker, **kwargs):

    # Outline
    # if composer issubclass(checker):
    # elif composer issubclass(dict):
    # else (none):

    if composer is not None:
        # print(type(composer))
        print('composer isinstance dict')
        if not issubclass(type(composer), checker):
            print('composer is not instance checker')
            if project_path:
                instance.project_path = Path(project_path) / Path(instance.project)
            else:
                instance.project_path = Path(os.getcwd()) / Path(instance.project)
            Path.mkdir(instance.project_path, parents=True, exist_ok=True)
            print('1project_path=%s' % instance.project_path)
            instance.removed_pm_config(kwargs)
        else:
            setattr(composer, 'project', project)
            for key, value in composer.__dict__.items():
                setattr(instance, key, value)
                print('key:' + str(key) + '\nvalue: ' + str(value))
            if 'project_path' not in composer.__dict__.keys():
                if project_path:
                    instance.project_path = instance.repo_path
                else:
                    instance.project_path = Path(os.getcwd()) / Path(instance.project)
            print('2project_path=%s' % instance.project_path)
    return instance


def removed_pm_config(instance, project, project_path, basic_project=True, custom=None):
    # Parameters:
        # P1 Do we use a basic project cookie?
        # P2 Do we use custom file locations?
            # Use a dictionary with the names below
        # If we are using P1 and P2 then P2 will override P1.
            # Below are the variables that can be overridden
            # The key is one of those variable names
            # The value is a Path-Like object
    if basic_project:
        Kitchen = Oven(project=project, basic_project=basic_project)
        Kitchen.bake_the_project(cookie_jar=project_path)

    instance.project_index = project_path / Path('index')
    instance.user_db = project_path / Path('databases')
    instance.ncbi_db_repo = instance.user_db / Path('NCBI')
    instance.project_database = instance.user_db / Path(project)
    instance.raw_data = project_path / Path('raw_data')
    instance.data = project_path / Path('data')
    instance.research_path = project_path

    if custom:
        for key, value in custom.items():
            setattr(instance, key, value)

    Path.mkdir(instance.project_index, exist_ok=True)
    Path.mkdir(instance.user_db, exist_ok=True)
    Path.mkdir(instance.ncbi_db_repo, exist_ok=True)
    Path.mkdir(instance.project_database, exist_ok=True)
    Path.mkdir(instance.raw_data, exist_ok=True)
    Path.mkdir(instance.data, exist_ok=True)