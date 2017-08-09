# Author: ${AUTHOR}
# Date Created: {$DATE}
# Project Name: QC_test
# Description: Testing the quality control of the new phylogenetics pipeline.

#PBS -S /bin/bash
#PBS -m bea
#PBS -M ${EMAIL}
#PBS -l mem=${GBS}
#PBS -N ${AUTHOR}_${GENE}
#PBS -o ${AUTHOR}_${GENE}.o$PBS_JOBID
#PBS -e ${AUTHOR}_${GENE}.e$PBS_JOBID
#PBS -j oe

# change into ${RAWDATA}
cd ${RAWDATA}

python36 ${PYTHONFILE} ${GENE} ${RAWDATA}

mail -s "${AUTHOR}_${GENE} script completed" ${EMAIL} <<< 'Check your output'