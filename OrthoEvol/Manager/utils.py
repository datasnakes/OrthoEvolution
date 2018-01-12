import yaml


def parse_db_config_file(config_file):
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