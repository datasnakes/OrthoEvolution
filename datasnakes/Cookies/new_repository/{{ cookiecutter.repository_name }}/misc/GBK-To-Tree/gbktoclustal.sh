#!/usr/bin/csh
# -*- coding: utf-8 -*-
# This is a simple bash script to run the python script
alias d='date "+DATE: %m/%d/%y%nTIME: %r-%Z"'
alias logdate='date +%m-%d-%y@%I:%M-%p.log'
logname=$(logdate)
echo "#--------------------------------------------------------------" >> $logname.log
echo Script started at $(d) >> $logname.log
echo "#--------------------------------------------------------------" >> $logname.log
python GBK-To-Clustal.py >> $logname.log
echo "#--------------------------------------------------------------" >> $logname.log
echo Script finished at $(d) >> $logname.log
echo "#--------------------------------------------------------------" >> $logname.log
mail -s "MCSR LOG: GBK-To-Clustal Script Update" shutchins2@umc.edu < $logname.log
mv $logname.log /work5/r2295/bin/LogFiles/GBK-to-Clustal