# -*- coding: utf-8 -*-
from Datasnakes.Tools.sge import SGEJob

myjob = SGEJob(base_jobname='testjob', email_address='shutchins2@umc.edu')

code = "print('hi')"
myjob.submit(code)
