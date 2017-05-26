"""Test for flask integration."""
import shutil
import os
from datasnakes.Manager.utils.mana import Mana, WebMana
x = Mana(repo='Tester', new_repo=True)
y = WebMana(repo='Tester', website='Vall', new_website=True)

# Delete the created `Tester` directory
if os.path.isdir('Tester') == True:
    print('Tester directory exists. Now delete it.')
    shutil.rmtree('Tester')
