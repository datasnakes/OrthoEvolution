# -*- coding: utf-8 -*-
from baseftp import BaseFTP
from ftplib import FTP
from Datasnakes import DatasnakesDevelopmentWarning


class NCBIFTP(BaseFTP):
    """Access NCBI's FTP servers with ease."""
    # raise DatasnakesDevelopmentWarning('NCBIFTP is not ready to be used.')
    def __init__(self, email, ftpsite='ftp.ncbi.nlm.nih.gov'):
        BaseFTP.__init__(self, email, ftpsite)


# TODO Create function to download blast databases
# TODO Create function to list databases
