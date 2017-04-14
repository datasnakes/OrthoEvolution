
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
PACKAGES = [
    'lib',
    'lib.scripts',
    'lib.scripts.biosql',
    'lib.scripts.blast',
    'lib.scripts.ftp',
    'lib.scripts.genbank',
    'lib.scripts.manager',
    'lib.scripts.multiprocessing',
    'lib.scripts.phylogenetic_analyses'
]
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='orthologs',
    description="A project that will help to analyze orthologous gense.",
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
    packages=
)
