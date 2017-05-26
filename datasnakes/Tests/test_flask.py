"""Test for flask integration."""
import shutil
import os
from datasnakes.Manager import Mana, WebMana

def flask_test(directoryname='Test-Directory', website='TestSite'):
    Mana(repo=directoryname, new_repo=True)
    WebMana(repo=directoryname, website=website, new_website=True)

    # Delete the created `Tester` directory
    if os.path.isdir(directoryname) == True:
        print('%s directory exists. Now delete it. Test passed.' % directoryname)
        shutil.rmtree(directoryname)
    else:
        print('%s directory DOES NOT exists. Debug failed test.' % directoryname)