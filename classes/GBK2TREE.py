# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 15:47:16 2017

@author: shutchins2
"""

# Classes and/or definitions to be used for GBK Feature Extraction,
# alignment, PhyML, PAML

import csv


class CreateLists:
    """ Use this class to create various lists need for the related scripts."""

    def org_list(file):

        olist = []  # Initialize list of organisms
        olist.append('')
        o = open(file)  # Open a comma delimited list of organisms.
        x = csv.reader(o)
        for org in x:  # Format a list of organisms
            org = str(org)
            org = org.replace("'", "")
            org = org.replace("[", "")
            org = org.replace("]", "")
            org = org.replace(" ", "_")
            olist.append(org)
        print("This is the list of organisms: " + "\n")
        print(olist)  # Print the list of organisms

        # To do - [] Complete these definitions