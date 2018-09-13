from pathlib import Path
from OrthoEvol.Cookies.cookie_jar import Oven
from OrthoEvol.Tools.logit import LogIt

class OrthologUtils(BlastUtils, GenbankUtils):

    def __init__(self):
        BlastUtils.__init__(self)
        GenbankUtils.__init__(self)

    def attribute_config(self, cls, composer, checker, project=None, project_path=None, checker2=None):
        """Set/configure attributes.

        Attribute Configuration takes an instance of a class and sets various
        attributes. The attributes are set by determining the type of configuration.
        The class attributes can be composed by another class, they can be set with
        a dictionary, or they can be set using the basic project template.

        :param cls: An instance of a class that will retain the attributes.
        :param composer: A class that will yield attributes to the cls parameter.
        :param checker: A checker class used to check the type of the composer.
                Dictionary composers will be treated differently.
        :param project:  The name of the project.
        :param project_path:  The relative path of the project.
        :return:  Returns the instance (cls) with new attributes.
        """
        ac_log = LogIt().default(logname="%s" % cls.__class__.__name__, logfile=None)
        if checker2:
            check2 = issubclass(type(composer), checker2)
        else:
            check2 = None
        # Attribute configuration using checker composition.
        if issubclass(type(composer), checker) or check2:
            for key, value in composer.__dict__.items():
                setattr(cls, key, value)
            ac_log.info("The attribute configuration was accomplished by composing %s with %s." % (cls.__class__.__name__, composer.__class__.__name__))

        # Attribute configuration using a dictionary.
        elif isinstance(composer, dict):
            for key, value in composer.items():
                setattr(cls, key, value)
            ac_log.info("The attribute configuration of %s was accomplished by using a dictionary." % cls.__class__.__name__)

        # Attribute configuration without composer
        elif composer is None:
            if not (project or project_path):
                raise BrokenPipeError("Without the Project Management class, a project name and project path must be included.")
            cls = self.standalone_config(cls, project, project_path)
            ac_log.info("The attribute configuration of %s was accomplished by using a standalone project." % cls.__class__.__name__)
        # Make sure self.project and self.project_path have values
        if not (cls.project or cls.project_path):
            raise BrokenPipeError("The project name and project path attributes have not been set.")

        return cls

    def standalone_config(self, cls, project, project_path, new=False, custom=None):
        """Configure a standalone project.

        A standalone configuration uses the variables listed.  These variables are
        either mapped to a basic project, used in a custom configuration, or they
        are mapped to some basic project directories with some custom options.

        :param cls: An instance of a class that will retain the attributes.
        :param project: The name of the project.
        :param project_path: The relative path of a project.
        :param new: The new project flag.
        :param custom: The custom flag which can be None or a dictionary.
        :return: Returns the instance (cls) with new attributes.
        """

        cls.project = project
        cls.project_path = project_path / Path(project)
        cls.project_index = cls.project_path / Path('index')
        cls.user_index = cls.project_path / Path('index')
        cls.db_archives = cls.project_path / Path('archive')
        cls.raw_data = cls.project_path / Path('raw_data')
        cls.data = cls.project_path / Path('data')
        cls.research_path = cls.project_path
        cls.user_db = cls.project_path / Path('databases')
        cls.project_database = cls.user_db / Path(project)
        cls.itis_db_repo = cls.user_db / Path('ITIS')
        cls.ncbi_db_repo = cls.user_db / Path('NCBI')
        cls.blast_db = cls.ncbi_db_repo / Path('blast') / Path('db')
        cls.windowmaker_files = cls.ncbi_db_repo / Path('blast') / Path('windowmaker_files')
        cls.ncbi_taxonomy = cls.ncbi_db_repo / Path('pub') / Path('taxonomy')
        cls.NCBI_refseq_release = cls.ncbi_db_repo / Path('refseq') / Path('release')

        # Use the basic_project cookie to create the directory structure
        if new or (not Path(cls.project_path).is_dir()):
            Kitchen = Oven(project=project, basic_project=True)
            Kitchen.bake_the_project(cookie_jar=project_path)

        # Use the custom dictionary to set the path variables
        # and to make the directories if necessary.  This overrides
        if custom:
            for key, value in custom.items():
                setattr(cls, key, value)
                if not Path(str(value)).is_dir():
                    Path.mkdir(value, exist_ok=True)

        return cls
