#!/bin/sh
# -*- coding: utf-8 -*-
# This is a simple bash script that performs a NCBI Blast to GBK Feature extraction.

alias d='date "+DATE: %m/%d/%y%nTIME: %r-%Z"'
alias logdate='date +%m-%d-%y@%I:%M-%p.log'
logname=$(logdate)

echo "#--------------------------------------------------------------" >> blast2fasta_$logname
echo Script started at $(d) >> blast2fasta_$logname
echo "#--------------------------------------------------------------" >> blast2fasta_$logname

echo "Run GetGiLists.py" >> blast2fasta_$logname
python3 GetGiLists.py >> blast2fasta_$logname
echo "#--------------------------------------------------------------" >> blast2fasta_$logname

echo "Run AccCollect.py" >> blast2fasta_$logname
python3 AccCollect.py >> blast2fasta_$logname
echo "#--------------------------------------------------------------" >> blast2fasta_$logname

echo "Run ScriptName.py" >> blast2fasta_$logname
python3 script_name.py >> blast2fasta_$logname
echo "#--------------------------------------------------------------" >> blast2fasta_$logname

echo Script finished at $(d) >> blast2fasta_$logname
echo "#--------------------------------------------------------------" >> blast2fasta_$logname

mail -s "Blast2Fasta Update" shutchins2@umc.edu < blast2fasta_$logname

mv $logname /work2/vallender/Projects/GPCR-Orthologs/DocsAndFiles/LogFiles

# add Rob & Dr. Vallender to the mail output using -c rgilmore@umc.edu -c evallender@umc.edu
# >> appends to the text file
# < file.txt sends the text file to the corresponding email addresses
## <<< 'hey' sends text to the mail command as the body