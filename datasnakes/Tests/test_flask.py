"""Test for flask integration."""
import shutil
import os
from datasnakes.Manager import Mana, WebMana


def flask_test(directoryname='Test-Directory', website='TestSite'):
    """"Test flask by creating a website/directory."""
    Mana(repo=directoryname, new_repo=True)
    WebMana(repo=directoryname, website=website, new_website=True)
    # TODO-SDH Write a system exit for the flask tests.

    # Delete the created `Tester` directory
    if os.path.isdir(directoryname) is True:
        print('%s directory exists. Delete it. Test passed.' % directoryname)
        shutil.rmtree(directoryname)
    else:
        print('%s directory DOES NOT exists. Test failed.' % directoryname)
