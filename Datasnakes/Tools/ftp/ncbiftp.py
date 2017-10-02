# -*- coding: utf-8 -*-
from .baseftp import BaseFTP
from ftplib import FTP
from Datasnakes import DatasnakesDevelopmentWarning


class NCBIFTP(BaseFTP):
    """Access NCBI's FTP servers with ease."""
    raise DatasnakesDevelopmentWarning('NCBIFTP is not ready to be used.')


# TODO Create function to download blast databases
# TODO Create function to list databases
