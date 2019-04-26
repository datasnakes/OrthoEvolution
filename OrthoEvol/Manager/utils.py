import yaml
from OrthoEvol.Tools.sge import SGEJob


class ManagerUtils(object):
    def __init__(self):
        super().__init__()

    def parse_db_config_file(self, config_file):
        kw ={}
        db_config_strategy = {}
        with open(config_file, 'r') as cf:
            db_config = yaml.load(cf)
            # Get the configuration for the desired strategy
            for key, value in db_config["Database_config"].items():
                if isinstance(value, dict):
                    db_config_strategy[key] = value
                    continue
                # Get the parameters for the Base class
                else:
                    kw[key] = value
        return db_config_strategy, kw

    def refseq_jobber(self, email_address, base_jobname, id, code):
        job = SGEJob(email_address=email_address, base_jobname=base_jobname % str(id))
        job.submit_pycode(code=code, wait=False, cleanup=False)

    # def template_jobber(email_address, base_jobname, id, code):
    #     job = SGEJob(email_address=email_address, base_jobname=base_jobname)
    #     job.submit_pycode(code=code, wait=True, cleanup=True)
