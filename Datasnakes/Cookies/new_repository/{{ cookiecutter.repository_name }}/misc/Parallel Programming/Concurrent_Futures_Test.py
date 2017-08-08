# -*- coding: utf-8 -*-
"""
Date created: Thu Mar  2 12:03:42 2017
Author: S. Hutchins

Script description:

"""
import time
try:
    from urllib.request import urlopen
except:
    from urllib2 import urlopen
from multiprocessing.dummy import Pool as ThreadPool

urls = ['https://www.google.co.uk', 'https://facebook.com', 'http://news.bbc.co.uk', 'http://www.reddit.com']

start = time.time()
pool = ThreadPool(4)
results = pool.map(urlopen, urls)
pool.close()
pool.join()
print ("Elapsed time: %s" % (time.time()-start))