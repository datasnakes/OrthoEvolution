from Datasnakes.Manager.utils.mana import ProjMana
from Datasnakes.Orthologs.Blast.blastn import BLASTn
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Align.alignment import Alignment

project='test-project'

pm = ProjMana(repo='test-repo', user='test-user', project=project, research='test-research', research_type='comparative_genetics',
              new_repo=True, new_user=True, new_project=True, new_research=True)

bn = BLASTn(project=project, copy_from_package=True, MAF='MAFV3.2.csv', proj_mana=pm)

gb = GenBank(project=project, blast=bn)
al = Alignment(project=project, aln_program='GUIDANCE2', genbank=gb)
al.align(seqFile='test.ffn', msaProgram='CLUSTALW', seqType='nuc')
