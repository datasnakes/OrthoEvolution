"""Test for flask integration."""
import shutil
import os
from Datasnakes.Manager.utils import Mana, WebMana


def flask_test(directoryname='Flask-Test', website='TestSite'):
    """"Test flask by creating a website/directory."""
    try:
        Mana(repo=directoryname, new_repo=True)
        WebMana(repo=directoryname, website=website, new_website=False)
        # new_website=False in order to prevent starting of server.

        if os.path.isdir(directoryname) is True:
            print('%s directory exists. Deleted it. Test passed.' % directoryname)
            # Delete the created `Tester` directory
            shutil.rmtree(directoryname)
        else:
            print('%s directory was not created.' % directoryname)
    except:  # Raise an exception if the directories are not created.
        raise Exception('Error with Mana or WebMana.')
