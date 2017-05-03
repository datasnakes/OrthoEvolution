##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
Orthologs-Project
proj_mana.py updated on 4/3/2017 at 3:39 PM
##############################################################################

    Input:

    Output:

    Description:

##############################################################################
@author: rgilmore
"""
##############################################################################
# Libraries:

import os
#from dir_mana import dir_mana
#from lister import Lister
import json
import tablib
##############################################################################
# Custom Class Initializations
# :
# Use directory_management() class here so that we can stay organized
# and more easily access the proper directories on command
home = os.getcwd()
project = "Orthologs-Project"
user = "rgilmore"
# Use lister() class here so that we can easily access our Master RNA Accession File
#what = Lister('MAFv3.1.csv')  # Always make sure this file name is correct

## Add a path that contains custom libraries for import
#os.sys.path.append()
##############################################################################
# Global Initializations:

##############################################################################

class proj_mana(object):

    def __init__(self, projects='', datasets='', research_targets='', js_load=False):
        data_config = tablib.Dataset().load(open('data_config.yaml').read())

        self.projects = projects
        self.datasets = datasets
        self.research_targets = research_targets
        if js_load is True:
            self.proj_mana_json = {}
            with open('proj_mana.json', 'r') as file:
                self.proj_mana_json = json.load(file)[0]
        self.manage(self.projects, self.datasets, self.research_targets, js_load)


    def yes_or_no(self, question):
        reply = str(input(question + ' (y/n): ')).lower().strip()
        if reply[0] == 'y':
            return True
        if reply[0] == 'n':
            return False
        else:
            return self.yes_or_no("Invalid answer... please enter ->")

    def ui(self, handle_key, handle_value):
        print('Here are the current %s\'s that we have:' % handle_key)
        for items in handle_value:
            print(items)
        adding = self.yes_or_no("\nWould you like to add more %s\'s?" % handle_key)
        names = []
        while adding is True:
            name = input("\nAdding a new %s....\nWhat's the name of the %s?\n" % (handle_key, handle_key))
            valid = self.yes_or_no("%s? Is that correct?" % name)
            if valid is True:
                names.append(name)
            else:
                print('%s was not added to the %s list.' % (name, handle_key))
                continue
            ans = self.yes_or_no("Would you like to add more data-sets?")
            if ans is True:
                adding = True
            else:
                adding = False

    def manage(self, projects, datasets, r_targets, js_load):
        # //TODO-ROB Make a JSON dump file for the projects/datasets/targets dictionary
        # //TODO-ROB Make a custom SQL database with project/directory management information
        if js_load is False:
            # Initialize the project information
            projects = list(projects) + ['Orthologs-Project', 'RNAseq-Project']
            projects = list(set(projects))
            r_targets = list(r_targets) + ['comparative_evolution', 'comparative_polymorphism', 'natural_selection']
            r_targets = list(set(r_targets))
            datasets = list(datasets) + ['GPCR_dataset', 'KARG_dataset', 'Hall_dataset']
            datasets = list(set(datasets))

        # Ask the user for new project information
        # Iterate through project data asking for user input along the way and return a dictionary
        ds_rt = {}  # dict {data-set : research_targets}
        p_ds = {}  # dict {project : ds_rt}
        print("Choose a research target for each project")
        for proj in projects:
            print('\nProject: %s' % proj)
            for ds in datasets:
                print('Data-set: %s' % ds)
                valid = self.yes_or_no("Is the %s data-set being used in the %s project?" % (ds, proj))
                if valid is False:
                    continue
                targs = []
                for target in r_targets:
                    print('Research Target: %s' % target)
                    valid = self.yes_or_no('Is the %s research target being used to analyze the %s dataset?' % (target, ds))
                    if valid is True:
                        targs.append(target)
                    else:
                        continue
                ds_rt[ds] = targs
            p_ds[proj] = ds_rt
        with open("proj_mana.json", 'w') as file:
            json.dump(cl_p, file)
        return p_ds

