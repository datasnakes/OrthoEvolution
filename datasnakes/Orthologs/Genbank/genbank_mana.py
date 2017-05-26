import os
from pathlib import Path
from datasnakes.Orthologs.Genbank.genbank import GenBank as GB
from datasnakes.Tools.ftp import ftp2db as FTP


class GenBankMana(FTP, GB):

    def __init__(self):
        FTP.__init__()
        GB.__init__()
        print("I don't do anything but look cool.")
