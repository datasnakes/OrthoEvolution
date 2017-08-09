#import os
#from pathlib import Path
from .genbank import GenBank as GB
from Datasnakes.Tools.ftp.ftp2db import FTP

class GenBankMana(FTP, GB):
    """
    """
    def __init__(self):
        FTP.__init__(self)
        GB.__init__(self)
        print("I don't do anything but look cool.")
