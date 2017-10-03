# -*- coding: utf-8 -*-
from baseftp import BaseFTPClient
from ftplib import FTP
from Datasnakes import DatasnakesDevelopmentWarning


class NCBIFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""
    # raise DatasnakesDevelopmentWarning('NCBIFTP is not ready to be used.')
    def __init__(self, email, ftpsite='ftp.ncbi.nlm.nih.gov'):
        super().__init__(ftpsite, email)


# TODO Create functions to download blast databases
# TODO Create function to list databases
