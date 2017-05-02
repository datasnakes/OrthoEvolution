from manager.lister import Lister
from manager.BLASTingTemplate import BLASTingTemplate
from manager.BlastnScript.BLASTn import BLASTn

x = BLASTn('MAFV3.2.csv')  # This is a template for GPCR project

BLASTER = x.blast_config

BLASTER(x.blast_human, 'Homo sapiens', auto_start=True)
