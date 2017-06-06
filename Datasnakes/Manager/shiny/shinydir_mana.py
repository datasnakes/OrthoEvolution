"""Directory management for shiny app."""
import os
import platform
from dir_mana import dir_mana

where = dir_mana(os.path.dirname(os.path.abspath(__file__)), 'GPCR-Orthologs-Project')
os.chdir(where.APP_DATA)
computer = str(platform.uname()[1])

if computer == 'gpcr-orthologs-project':
    with open('shinydir_mana.txt', 'w') as file:
        for key, value in where.__dict__.items():
            file.write('%s#%s\\\n' % (key, value))
elif 'UMC' in computer:
    with open('shiny_dir_mana.txt', 'w') as file:
        for key, value in where.__dict__.items():
            file.write('%s%s\\\n' % (key, value))
