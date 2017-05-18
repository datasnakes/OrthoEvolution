from Orthologs.blast.blastn import BLASTn
from Manager.utils import mana
import shutil
from pathlib import Path
# Initializations
repo = "Test1"
user = "rgilmore"
project = "BLASTest"
research = "comparative_genetics"
research_type = "public"
# Directory setup now combined with blast setup
x = mana.ProjMana(repo, user, project, research, research_type, new_project=True, new_research=True, new_repo=True, new_user=True)
shutil.copy(str(x.index / Path('MAFV3.2.csv')), str(x.project_index))
#x = BLASTn(repo, user, project, research, research_type, template="MAFV3.2.csv", new_project=True, new_research=True, new_repo=True, new_user=True)
x = BLASTn(repo, user, project, research, research_type, template="MAFV3.2.csv")
BLASTER = x.blast_config
BLASTER(x.blast_human, 'Homo_sapiens', auto_start=True)
