"""BLASTn test"""
import shutil
#from pathlib import Path
from Datasnakes.Manager.utils.mana import ProjectManagement
from Datasnakes.Manager import index
from Datasnakes.Orthologs.Blast import BLASTn
# from ete3 import NCBITaxa
# NCBITaxa().update_taxonomy_database()
import pkg_resources
# Initializations
repo = "Test2"
user = "johndoe"
project = "BLASTest"
research = "comparative_genetics"
research_type = "public"
template = pkg_resources.resource_filename(index.__name__, "MAFV3.2.csv")
# Directory setup now combined with blast setup
# a = ProjectManagement(
#     repo=repo,
#     user=user,
#     project=project,
#     research=research,
#     research_type=research_type, new_repo=True)


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
    template=template, new_repo=True, new_project=True, new_research=True, new_user=True, copy_from_package=True, MAF='MAFV3.2.csv')
BLASTER = x.blast_config
BLASTER(x.blast_human, 'Homo_sapiens', auto_start=True)
# x.post_blast_analysis("Vall")
