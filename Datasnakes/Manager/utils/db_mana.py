from Datasnakes.Tools.ftp.ftp2db import FTP
from Datasnakes.Manager.utils.mana import UserMana as UM


class DbMana(FTP, UM):

    def __init__(self):
        FTP.__init__(self)
        UM.__init__(self)