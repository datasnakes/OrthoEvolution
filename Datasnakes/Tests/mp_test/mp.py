import subprocess
import re
print("Begin subprocess.run call")
qsub = subprocess.check_output(['qsub mp.sh'], universal_newlines=True)
job_id = re.findall(r'\d+', qsub)[0]

