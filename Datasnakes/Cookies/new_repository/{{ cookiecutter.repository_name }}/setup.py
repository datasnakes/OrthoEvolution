# Used:
# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://github.com/biopython/biopython/blob/master/setup.py
# TODO-ROB: 1.  Generate a Template BioSQL database using sqlite3
from setuptools import setup
# To use a consistent encoding
from codecs import open
from os import path

home = path.abspath(path.dirname(__file__))
name = 'Orthologs-Project'



# Get the long description from the README file
with open(path.join(home, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Set up an initial configuration file




setup(
    name=name,
    description="A project that will help to analyze orthologous genes.",
    version='0.1.0',
    long_description=long_description,
    url='https://github.com/robear22890/Orthologs-Project',
    license='?',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
        'Natural Language :: English'
    ],
    #packages=PACKAGES,
    install_requires=[],
)

