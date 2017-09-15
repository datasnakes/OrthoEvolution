# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 10:55:30 2016

@author: Shaurita D. Hutchins
"""
# This program will open a ssh terminal and login to the MCSR.

import os
import subprocess

os.chdir('C://Users//shutchins2//Desktop//PUTTY')
cmd = "plink -v -load MCSR -pw Medschool17"
process = subprocess.Popen(cmd)
stdout, stderr = process.communicate()
