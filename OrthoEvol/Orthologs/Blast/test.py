from OrthoEvol.Orthologs.Blast import OrthoBlastN, BaseBlastN

# This is more pythonic with YAML loading
blastconfig = {
    "project": "sdh-test",
    "blast_method": 3,
    "taxon_file": None,
    "go_list": None,
    "post_blast": True,
    "template": None,
    "save_data": True,
    "copy_from_package": True,
    "MAF": 'test_blast.csv',
    "project_path": None,
    "proj_mana": None
}


myblast = BaseBlastN(**blastconfig)
myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)
