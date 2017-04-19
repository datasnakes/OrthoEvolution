import tablib
import os
from pathlib import Path
import configparser


# Configure "environment variables"
class config_env(object):

    def __init__(self):
        count = Path.cwd().parts.__len__()
        i = 0
        while i < count:
            try:
                config = tablib.Dataset().load(open('init.yaml').read())
                break
            except FileNotFoundError:
                i += 1
                os.chdir('..')

        self.home = config['location']
        self.package = config['package']
