# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:43:50 2017

@author: shutchins2

install 'progressbar33' via pip
"""

import sys
import time

from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer

pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=300).start()
for i in range(300):
    time.sleep(0.01)
    pbar.update(i+1)
pbar.finish()