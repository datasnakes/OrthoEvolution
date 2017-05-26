# -*- coding: utf-8 -*-
"""
Date created: Thu Mar  2 12:22:09 2017
Author: S. Hutchins

Script description: Multiprocessing test for the MCSR.

"""
from multiprocessing import Pool

def count_chars(filename):
    with open(filename) as f:
        return len(f.read())

def success_callback(result):
    print('Get result', result)

def error_callback(e):
    print('Throw exception', e)

filelist = ['test1.txt', 'test2.txt', 'test3.txt']

if __name__ == '__main__':
    with Pool(processes=2) as pool:
        async_result = pool.map_async(count_chars, filelist,
                                     callback=success_callback,
                                     error_callback=error_callback)
        # Wait for the task to complete.
        async_result.wait()

