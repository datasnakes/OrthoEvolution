Naming Conventions and Paths *by class*
=======================================

CompGenAnalysis(ProjectMana) *comp\_gen.py*
-------------------------------------------

User Generated Files
~~~~~~~~~~~~~~~~~~~~

-  taxon\_file

   -  path = *project\_index* / Path(taxon\_file)
   -  name = "*projectname*\ \_tax.txt"

-  paml\_file

   -  path = *project\_index* / Path(paml\_file)
   -  name = "*projectname*\ \_paml.txt"
   -  TODO-ROB: Deprecate this file

-  acc\_file

   -  path *project\_index* / Path(acc\_file)
   -  name = "*projectname*\ \_acc.csv"

Programmatically Generated Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  log\_files

   -  1
   -  2

-  my\_gene\_file

   -  path = \_project\_data / Path(my\_gene\_file)
   -  name = "*projectname*\ \_mygene.csv"

-  post\_blast\_analysis\_file

   -  path *project\_index* / Path(pba\_file)
   -  name = "*projectname*\ \_pba.xlsx"

BLASTAnalysis(CompGenAnalysis) *ncbi\_blast.py*
-----------------------------------------------

User Generated Files
~~~~~~~~~~~~~~~~~~~~

-  template\_file

   -  path = *project\_index* / Path(template\_file)
   -  name = "*projectname*\ \_temp.csv"

Programmatically Generated Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  log\_files

   -  TODO-ROB Add loggers
   -  accession/time logger
   -  post blast analysis logger

-  building\_file

   -  path = *project\_raw\_data* / Path(builder\_file)
   -  name = "*projectname*\ \_building.csv"

-  building\_time\_file

   -  path = *project\_raw\_data* / Path(builder\_time\_file)
   -  name = "*projectname*\ \_building\_time.csv"

BLASTn(BT) *blastn.py*
----------------------

User Generated Files
~~~~~~~~~~~~~~~~~~~~

-  template file

   -  path = *project\_index* / Path(template file)
   -  name = "*users\_file\_name*.csv"

Programmatically Generated Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  gi list binary files

   -  path = *project\_index* / Path(binary\_files)
   -  names = "*TAXID*\ \_gi"

-  blast xml files

   -  path = *project\_raw\_data* / Path(gene) / Path(blast xml files)
   -  name = "*gene*\ \_\ *organism*.xml

-  master accession file

   -  path = *project\_data* / Path(MAF)
   -  name = "*projectname*\ \_MAF.csv"
   -  ***Note***: *Only exists if the BLAST finishes all the genes.*

Naming Conventions and Paths *by location*
==========================================
