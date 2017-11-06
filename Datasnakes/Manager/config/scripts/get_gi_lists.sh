# TODO: GI deprecation.
# #PBS -S /bin/bash
# #PBS -m bea
# #PBS -l select=10:ncpus=1:mem=2gb -l place=free
# #PBS -l cput=24:00:00
# #PBS -l walltime=32:00:00
# #PBS -N getgilists
# #PBS -o getgilists.o${PBS_JOBID}
# #PBS -e getgilists.e${PBS_JOBID}
# #PBS -j oe
#
# cd ${PBS_O_WORKDIR}
# mpiexec python get_gi_lists.py