from Orthologs.manager.blast.blastn import BLASTn

x = BLASTn('MAFV3.2.csv')  # This is a template for GPCR project

BLASTER = x.blast_config

BLASTER(x.blast_human, 'Homo_sapiens', auto_start=True)
