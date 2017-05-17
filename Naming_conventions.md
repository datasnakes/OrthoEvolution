#  Naming Conventions and File Placements

## CompGenAnalysis(ProjectMana) _comp_gen.py_
### User Generated Files
* taxon_file
    * path = _user_index_ / Path(taxon_file)
    * name = "_projectname_\_tax.txt"
* paml_file
    * path = _user_index_ / Path(paml_file)
    * name = "_projectname_\_paml.txt"
    * TODO-ROB:  Deprecate this file
* acc_file
    * path _user_index_ / Path(acc_file)
    * name = "_projectname_\_acc.csv"

### Programmatically Generated Files
* my_gene_file
    * path = _project_data / Path(my_gene_file)
    * name = "_projectname_\_mygene.csv"
* post_blast_analysis_file
    * path _project_data_ / Path(pba_file)
    * name = "_projectname_\_pba.xlsx"

## BLASTAnalysis(CompGenAnalysis)  _ncbi_blast.py_
### User Generated Files