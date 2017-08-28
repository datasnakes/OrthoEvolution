# Author: $AUTHOR
# Date Created: $DATE
# Project Name: QC_test
# Description: Testing the quality control of the new phylogenetics pipeline.

. ~/.bashrc
#PBS -S /bin/bash
#PBS -m bea
#PBS -M $EMAIL
#PBS -l mem=$GBS
#PBS -N $AUTHOR_$GENE
#PBS -o $AUTHOR_$GENE.o$PBS_JOBID
#PBS -e $AUTHOR_$GENE.e$PBS_JOBID
#PBS -j oe

# echo environment variables for output file
echo "Author is $AUTHOR"
echo "Date and time is $DATE"
echo "Email is $EMAIL"
echo "Memory allocation is $GBS"
echo "Target gene is $GENE"
echo "raw_data path is $RAWDATA"
echo "Python worker file is $PYTHONFILE"

# change into $RAWDATA
cd $RAWDATA

python36 $PYTHONFILE $GENE $RAWDATA

mail -s "$AUTHOR_$GENE script completed" $EMAIL <<< 'Check your output'