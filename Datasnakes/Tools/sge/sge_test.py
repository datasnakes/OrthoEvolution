# -*- coding: utf-8 -*-
from Datasnakes.Tools.sge import SGEJob

myjob = SGEJob(email_address='shutchins2@umc.edu')

code = "test.py"
myjob.submit(code)
