from Datasnakes.Tools.ftp import FTP2DB
from Datasnakes.Manager.utils import UserMana as UM


class DbMana(FTP2DB, UM):

    def __init__(self):
        FTP2DB.__init__(self)
        UM.__init__(self)
