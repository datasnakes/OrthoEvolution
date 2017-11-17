from OrthoEvol.Tools.sge import SGEJob

myjob = SGEJob(email_address='shutchins2@umc.edu')

code = "test.py"
myjob.submit_pycode(code)
