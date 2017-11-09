from Datasnakes.Orthologs.Blast import OrthoBlastN
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
myblast = OrthoBlastN(proj_mana=None, project="sdh-test", project_path=path, **Blast_config)
myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)
