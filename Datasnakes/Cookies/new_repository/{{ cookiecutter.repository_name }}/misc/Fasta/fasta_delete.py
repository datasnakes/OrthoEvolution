# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:33:06 2016

@author: rgilmore
"""

import fnmatch
import os

os.chdir('/work5/r2294/bin/NCBI_data/Delete Me/')
home = os.getcwd()
for file in os.listdir():
    if os.path.isdir('/work5/r2294/bin/NCBI_data/Delete Me/%s' % file) == False:
        continue
    if os.path.islink('/work5/r2294/bin/NCBI_data/Delete Me/%s' % file) == True:
        continue
    os.chdir('/work5/r2294/bin/NCBI_data/Delete Me/%s' % file)
    print(os.getcwd())
    for FILE in os.listdir():
        if fnmatch.fnmatch(FILE, '*.f*') == True:
            os.remove(FILE)
            print(FILE, 'File Deleted!')
    if os.path.isdir('/work5/r2294/bin/NCBI_data/Delete Me/%s' % file) == True:
        os.removedirs('/work5/r2294/bin/NCBI_data/Delete Me/%s' % file)
        print(file, 'Directory Deleted')
        continue
    os.chdir(home)