import getpass
from datetime import datetime as d
from collections import OrderedDict

class BasePBSJob(object):
    """Base class for PBS jobs."""

    def __init__(self, author=getpass.getuser(), project_name="OrthoEvol", description="This is a basic pbs job",
                 date_format='%a %b %d %I:%M:%S %p %Y', chunk_resources=None, cput='72:00:00', walltime='48:00:00',
                 job_name=None, pbs_work_dir=None, script_cmd=None, email=None, directive_list=None):

        # Set up commented script header
        self.author = author
        self.project_name = project_name
        self.description = description
        self.date_format = date_format

        # Set up PBS directives/attributes
        # resources
        if chunk_resources is None:
            chunk_resources = OrderedDict({
                "select": 1,
                "ncpus": 1,
                "memgb": "6gb"
            })

        resource_list = []
        for k, v in chunk_resources.items():
            if v is not None:
                resource_list.append("%s=%s" % (k, v))
        self.resource_str = ":".join(resource_list)
        self.cputime = cput
        self.walltime = walltime
        # other attributes
        self.job_name = job_name
        self.directive_list = directive_list

        # Set up PBS variables
        self.pbs_work_dir = pbs_work_dir
        self.script_cmd = script_cmd
        self.email = email
