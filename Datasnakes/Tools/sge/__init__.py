import platform
from .sgeutils import randomid, basejobids, import_temp, writecodefile, file2str
from .qstat import Qstat
from .sgejob import SGEJob
from .sgeconfig import __DEFAULT__, __CUSTOM__
from .sgepipelinetask import SGEPipelineTask


if "windows" in platform.system().lower():
    raise ImportError("This module is only supported on linux/osx.")