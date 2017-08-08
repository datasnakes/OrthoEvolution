#PBS -S /bin/bash
#PBS -m bea
#PBS -M shutchins2@umc.edu
#PBS -l select=8:ncpus=1:mem=16gb -l place=free
#PBS -l cput=24:00:00
#PBS -l walltime=32:00:00
#PBS -N sdhtest
#PBS -o /work5/r2295/bin/MultiprocessingTest/sdhtest.o${PBS_JOBID}
#PBS -e /work5/r2295/bin/MultiprocessingTest/sdhtest.e${PBS_JOBID}
#PBS -j oe
cd ${PBS_O_WORKDIR}
#rm /work5/r2295/bin/MultiprocessingTest/sdhtest.o*
#rm /work5/r2295/bin/MultiprocessingTest/sdhtest.e*
python3 /work5/r2295/bin/MultiprocessingTest/mptest.py
echo "End of script."