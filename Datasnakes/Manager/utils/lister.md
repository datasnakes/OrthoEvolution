Lister
=======

Lister is a class that creates a Lister object, which grants access to all of the initial data that we will use for the project.

Lister uses pandas to manipulate our data so that we can more easily call it other processes.

## Description

The data to be initialized includes the following:

* Master_Accession_File.csv (include project name)
* common_names.csv (PAML) (rename)
* taxon_ids.txt

 **Add some examples of these files as template**

## Usage
Using docstrings for help:
```python
from manager.lister import Lister
help(Lister)
```
Using this class is simple:
```python
from lister import Lister
from pprint import pprint
x = Lister(acc_file = 'Master_Accession_File.csv', paml_file = 'common_names.csv', taxon_file = 'taxon_ids.txt')
```

### Parameters
```python
from lister import Lister
x = Lister(acc_file, paml_file, taxon_file, go_list = None, hgnc_file = False)
```
* **acc_file** (*'Master_Accession_file.csv'*) - Contains accession numbers for a group of genes ranked by tier.  Each gene has a group of
orthlogs used in our phylogenetic anlysis.
* **paml_file** (*'commonnames.csv'*) - Contains a list of shortened organism names used in the MSA files.  This is done to compy with
PAML.
* **taxon_file** (*'taxon_ids.txt'*) - Contains an ordered list of taxon ids
* **go_list** (*[[gene.1, org.1], ... , [gene.n, org.n]]*) - A nested list that can be used to get information about specific gene/org pairs.
* **hgnc_file** - For future implementation.  Used as a file handle to parse an HGNC *.csv* file.

#### Returned Variables
* x.gene_count
* x.org_count
* x.paml_org_list
* x.taxon_ids

### Lists
List that contain header info(**x.header**):
```python
pprint(x.header)
['Tier',
 'Gene',
 'Homo_sapiens',
 'Macaca_mulatta',
 'Mus_musculus',
 'Rattus_norvegicus',
...
 'Trichechus_manatus_latirostris',
 'Tupaia_chinensis',
 'Tursiops_truncatus']
```
List of Accessions (**x.acc_list**):
```python
pprint(x.acc_list)
['NM_000680.3',
 'NM_000679.3',
 'NM_000678.3',
...
 'xm_004368425.2',
 'XM_006155397.2',
 'XM_004317686.1']
```
List of Genes (**x.gene_list**):
```python
pprint(x.gene_list)
['ADRA1A',
 'ADRA1B',
 'CHRM2',
 ...
 'SSTR2',
 'TSHR',
 'VIPR1']
```
List of Organisms (**x.org_list**):
```python
pprint(x.org_list)
['Homo_sapiens',
 'Macaca_mulatta',
 'Mus_musculus',
 'Rattus_norveg',
...
'Trichechus_manatus_latirostris',
 'Tupaia_chinensis',
 'Tursiops_truncatus']
```
### Dictionaries
To use a dictionary of accessions.
```python
>>> pprint(x.acc_dict)
{'NM_000115.3': ['EDNRB', 'Homo_sapiens'],
 'NM_000145.3': ['FSHR', 'Homo_sapiens'],
 'NM_000164.3': ['GIPR', 'Homo_sapiens'],
...
 'NM_001001620.1': ['CCR3', 'Sus_scrofa'],
 'NM_001002911.3': ['GPR139', 'Homo_sapiens'],
 'NM_001002944.1': ['ADORA2B', 'Canis_lupus_familiaris']}
```
Dictionary of Genes is a nested dictionary. (**x.gene_dict**, **x.tier_dict**):
```python
pprint(x.gene_dict['HTR1A'])
{'Ailuropoda_melanoleuca': 'XM_002926305.1',
 'Bos_taurus': 'XM_600535.5',
 'Callithrix_jacchus': 'XM_008992005.2',
 ...
 'Tier': '1',
 'Trichechus_manatus_latirostris': 'xm_004374552.2',
 'Tupaia_chinensis': 'xm_006156821.1',
 'Tursiops_truncatus': 'xm_004325159.1'}

pprint(x.tier_dict['HTR1A'])
'1'
```
Dictionary of Organisms is a nested dictionary (**x.org_dict**)
```python
homosapiens_query = x.org_dict['Homo_sapiens'].values()
homosapiens_gene_list = x.org_dict['Homo_sapiens'].keys()
>>> pprint(list(homosapiens_query))
['NM_000680.3',
 'NM_000679.3',
 'NM_000678.3',
...
 'XM_011517263.2',
 'NM_000369.2',
 'NM_004624.3']
```

### Pandas Dataframes

Dataframe that uses Gene as an index (**x.df**):
```python
pprint(x.df.T.HTR1A)
['Tier                                            1
Homo_sapiens                          NM_000524.3
Macaca_mulatta                     NM_001198700.1
Mus_musculus                          NM_008308.4
...
Tupaia_chinensis                   xm_006156821.1
Tursiops_truncatus                 xm_004325159.1
Name: HTR1A, dtype: object']
```
Pivot Table MultiIndexed with pandas(**x.pt**):
```python
# #### Format the main pivot table #### #
self.pt = pd.pivot_table(pd.read_csv(self.__filename_path), index=['Tier', 'Gene'], aggfunc='first')
array = self.pt.axes[1].tolist() # Organism list
self.pt.columns = pd.Index(array, name='Organism')
```

If your data has tiers or is divided into groups (**x.get_tier_frame**, **x.tier_frame_dict**):
```python
Tiers = x.get_tier_frame('1')
Tiers.keys()
dict_keys(['1'])
Tiers = x.tier_frame_dict()
Tiers.keys()
dict_keys(['1', '2', '3', 'None'])
```

### Methods
Lookup Accessions (**x.get_accession(gene, org)**, **x.get_accesions(go_list=None)**)
```python
x.get_accession('HTR1A', 'Homo_sapiens')
'NM_000524.3'

go_list = [['HTR1A', 'Homo_sapiens'], ['HTR1A', 'Macaca_mulatta']]

x.get_accessions(go_list = go_list)
['NM_000524.3', 'NM_001198700.1']
```
Lookup a list of Accession for alignment(**x.get_accession_alignment(gene)**):
```python
pprint(x.get_accession_alignment('HTR1A'))
['NM_000524.3',
 'NM_001198700.1',
...
 'xm_006156821.1',
 'xm_004325159.1']
```
Get the master lists from a new dataframe(**self.get_master_list(df)**):
```python
from Manager.lister import Lister
import os
from pathlib import Path
y = Lister()
y.get_master_lists(csv_file='MAFV3.1.csv')
```