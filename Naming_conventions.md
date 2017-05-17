#  Naming Conventions and File Placements

## CompGenAnalysis(ProjectMana) _comp_gen.py_
### User Generated Files
* taxon_file
    * path = _project_index_ / Path(taxon_file)
    * name = "_projectname_\_tax.txt"
* paml_file
    * path = _project_index_ / Path(paml_file)
    * name = "_projectname_\_paml.txt"
    * TODO-ROB:  Deprecate this file
* acc_file
    * path _project_index_ / Path(acc_file)
    * name = "_projectname_\_acc.csv"

### Programmatically Generated Files
* log_files
    * 1
    * 2
* my_gene_file
    * path = _project_data / Path(my_gene_file)
    * name = "_projectname_\_mygene.csv"
* post_blast_analysis_file
    * path _project_index_ / Path(pba_file)
    * name = "_projectname_\_pba.xlsx"

## BLASTAnalysis(CompGenAnalysis)  _ncbi_blast.py_
### User Generated Files
* template_file
    * path = _project_index_ / Path(template_file)
    * name = "_projectname_\_temp.csv"

### Programmatically Generated Files
* log_files
    * TODO-ROB Add loggers
    * accession/time logger
    * post blast analysis logger
* building_file
    * path = _project_raw_data_ / Path(builder_file)
    * name = "_projectname_\_building.csv"
* building_time_file
    * path = _project_raw_data_ / Path(builder_time_file)
    * name = "_projectname_\_building_time.csv"

## BLASTn(BT) _blastn.py_
### User Generated Files
* template file
### Programmatically Generated Files
* gi list binary files
    * path = _project_index_ / Path(binary_files)
    * names = "_TAXID_\_gi"
*