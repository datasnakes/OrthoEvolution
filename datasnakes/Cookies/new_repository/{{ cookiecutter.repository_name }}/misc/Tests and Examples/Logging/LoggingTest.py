# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 18:12:54 2016

@author: shutchins2
"""

# Scipt that tests logging

import logging
FORMAT = '%(asctime)s | %(levelname)s: %(message)s'
logging.basicConfig(filename='example.log', format=FORMAT, datefmt='Date: %m/%d/%Y Time: %I:%M:%S %p', level=logging.INFO)
logging.info('Start the script.')
logging.warning('Stop here.')
