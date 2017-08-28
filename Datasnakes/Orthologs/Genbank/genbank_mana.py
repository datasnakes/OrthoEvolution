#import os
#from pathlib import Path
from .genbank import GenBank as GB

class GenBankMana(GB):
    """
    """
    def __init__(self):
        GB.__init__(self)
        print("I don't do anything but look cool.")
