from OrthoEvol.Orthologs.Blast import OrthoBlastN, BaseBlastN
import os

# This is more pythonic with YAML loading
Blast_config = {
  "taxon_file": None,
  "go_list": None,
  "post_blast": True,
  "template": None,
  "save_data": True,
  "copy_from_package": True,
  "MAF": 'test_blast.csv'
}


path = os.getcwd()
myblast = BaseBlastN(proj_mana=None, project="sdh-test", blast_method=None, project_path=path, **Blast_config)
myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)
