from .Cookies.utils import CookieUtils
from .Manager.utils import ManagerUtils
from .Orthologs.utils import OrthologUtils
from .Tools.otherutils.other_utils import OtherUtils


class FullUtilities(CookieUtils, ManagerUtils, OrthologUtils, OtherUtils):

    def __init__(self, ):
        CookieUtils.__init__(self)
        ManagerUtils.__init__(self)
        OrthologUtils.__init__(self)
        OtherUtils.__init__(self)
