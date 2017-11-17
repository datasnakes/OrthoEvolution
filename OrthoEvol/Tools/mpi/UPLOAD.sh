#!/bin/bash
#PBS -S /bin/bash
#PBS -m bea
#PBS -M rgilmore@umc.edu
#PBS -l select=8:ncpus=1:mem=16gb -l place=free
#PBS -l cput=24:00:00
#PBS -l walltime=32:00:00
#PBS -N robupload
#PBS -o /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/Lib/Log/robupload.o${PBS_JOBID}
#PBS -e /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/Lib/Log/robupload.e${PBS_JOBID}
#PBS -j oe
cd ${PBS_O_WORKDIR}
rm /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/robupload.o*
rm /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/robupload.e*
mpiexec python /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/multi_dbupload.py
echo "end"