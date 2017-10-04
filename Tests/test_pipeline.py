#from Datasnakes.Manager.utils import ProjectManagement
#from Datasnakes.Orthologs.Blast.blastn import CompGenBLASTn
#from Datasnakes.Orthologs.GenBank.genbank import GenBank
#from Datasnakes.Orthologs.Align.msa import MultipleSequenceAlignment
#
#project='test-project'
#
#pm = ProjectManagement(repo='test-repo', user='test-user', project=project,
#                       research='test-research', research_type='comparative_genetics',
#                       new_repo=True, new_user=True, new_project=True, new_research=True)
#
#bn = CompGenBLASTn(project=project, copy_from_package=True, MAF='MAFV3.2.csv', proj_mana=pm)
#
#gb = GenBank(project=project, blast=bn)
#al = MultipleSequenceAlignment(project=project, aln_program='GUIDANCE2', genbank=gb)
#al.align(seqFile='test.ffn', msaProgram='CLUSTALW', seqType='nuc')
#
## TODO Make test compatible for python3.4
