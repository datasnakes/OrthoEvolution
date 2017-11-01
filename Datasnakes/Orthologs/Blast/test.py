#from Datasnakes.Orthologs.Blast.utils import gi_list_config
#import os
#
#tids = ['9606', '9544', '9646', '37293', '9913', '9483', '9838', '9615']
#
#p = os.getcwd()
#
#gi_list_config(gi_list_path=p, taxonomy_ids=tids)

from Datasnakes.Orthologs.Blast import CompGenBLASTn
import os

# This is more pythonic with YAML loading
Blast_config = {
  "taxon_file": None,
  "go_list": None,
  "post_blast": True,
  "template": None,
  "save_data": True,
  "copy_from_package": True,
  "MAF": 'MAFV3.2.csv'
}


path = os.getcwd()
myblast = CompGenBLASTn(proj_mana=None, project="sdh-test", project_path=path, **Blast_config)
myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)
