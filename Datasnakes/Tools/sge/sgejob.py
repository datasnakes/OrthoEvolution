import contextlib
from subprocess import run, CalledProcessError, PIPE
import os

from Datasnakes.Tools.sge import basejobids, writecodefile, import_temp

class BaseJob(object):
    """Create a class for simple jobs."""
    def __init__(self):
        pass
    
    def _cleanup(base):
        os.remove(base + '.pbs')
        os.remove(base + '.py')


class SGEJob(BaseJob):
    """Create multiple jobs & scripts for each job to run based on
    splitting a list into chunks.
    """
    def __init__(self):
        pass
    
    def configure(self):
        pass
    
    def submit(self, code, default=True, prefix=None):
        """Creates and submit a qsub job. Also uses python code."""
        # TIP If python is in your environment as only 'python' update that.
        # TODO-SDH add a default parameter or custom parameter
        # If default, a python file will be created from code that is used.
        if self.default == default:
            baseid, base = basejobids()
            if prefix is not None:
                base = prefix + '_' + base
            writecodefile(filename=base, code=code, language='python')
            outfile = 'orthoevol_{}.out'.format(baseid)
            errfile = 'orthoevol_{}.err'.format(baseid)
            # Create the pbs script from the template or dict
            pbstemp = import_temp('temp.pbs')
        
            script_name = base.format(baseid)
        
            attributes = self.pbs_dict(outfile=outfile,
                                       errfile=errfile,
                                       script_name=script_name)
        
            with open(base + '.pbs', 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(attributes))
                pbsfile.close()
        else:
            raise NotImplementedError('Custom qsub jobs are forbidden.')
            # TODO Improve custom job creation
        #            pbstemp = self.import_temp('temp.pbs')
        #            with open(base + '.pbs', 'w') as pbsfile:
        #                pbsfile.write(pbstemp.substitute(self.pbs_dict))
        #                pbsfile.close()
        
        with contextlib.suppress(CalledProcessError):
            cmd = ['qsub ' + base + '.pbs']  # this is the command
            # Shell MUST be True
            cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        
            if cmd_status.returncode == 0:  # Command was successful.
                print('Job submitted.\n')
                # TODO add a check to for job errors or check for error file.
        
            else:  # Unsuccessful. Stdout will be '1'
                print("PBS job not submitted.")
#                _cleanup(base=base)
