from Orthologs.Blast.blastn import BLASTn
from Manager.utils import mana
import shutil
from pathlib import Path
#from ete3 import NCBITaxa
#NCBITaxa().update_taxonomy_database()

# Initializations
repo = "Test1"
user = "rgilmore"
project = "BLASTest"
research = "comparative_genetics"
research_type = "public"
# Directory setup now combined with blast setup
a = mana.ProjMana(repo=repo, user=user, project=project, research=research, research_type=research_type, new_project=True, new_research=True, new_repo=True, new_user=True)
shutil.copy(str(a.index / Path('MAFV3.2.csv')), str(a.project_index))
Path.mkdir(a.raw_data / Path('blast') / Path('gi_lists'), parents=True)
shutil.copy(str(a.index / Path('get_gi_lists.pbs')), str(a.raw_data / Path('blast') / Path('gi_lists')))
shutil.copy(str(a.index / Path('get_gi_lists.py')), str(a.raw_data / Path('blast') / Path('gi_lists')))
#x = BLASTn(repo, user, project, research, research_type, template="MAFV3.2.csv", new_project=True, new_research=True, new_repo=True, new_user=True)
x = BLASTn(repo=repo, user=user, project=project, research=research, research_type=research_type, template="MAFV3.2.csv")
BLASTER = x.blast_config
BLASTER(x.blast_human, 'Homo_sapiens', auto_start=True)
