"""BLASTn test"""
#import shutil
#from pathlib import Path
from Datasnakes.Manager import ProjMana
from Datasnakes.Orthologs.Blast import BLASTn

#from ete3 import NCBITaxa
# NCBITaxa().update_taxonomy_database()

# Initializations
repo = "Test2"
user = "johndoe"
project = "BLASTest"
research = "comparative_genetics"
research_type = "public"

# Directory setup now combined with blast setup
a = ProjMana(
    repo=repo,
    user=user,
    project=project,
    research=research,
    research_type=research_type)

# shutil.copy(str(a.index / Path('MAFV3.2.csv')), str(a.project_index))
# Path.mkdir(a.raw_data / Path('blast') / Path('gi_lists'), parents=True)
# shutil.copy(str(a.index / Path('get_gi_lists.sh')),
#             str(a.raw_data / Path('blast') / Path('gi_lists')))
# shutil.copy(str(a.index / Path('get_gi_lists.py')),
#             str(a.raw_data / Path('blast') / Path('gi_lists')))
x = BLASTn(
    repo=repo,
    user=user,
    project=project,
    research=research,
    research_type=research_type,
    template="MAFV3.2.csv")
BLASTER = x.blast_config
BLASTER(x.blast_human, 'Homo_sapiens', auto_start=True)
# x.post_blast_analysis("Vall")
