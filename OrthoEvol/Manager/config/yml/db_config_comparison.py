from pathlib import Path
# # webster config
#
# archive_options = {
#     "Full": Path(''),
#     "NCBI": Path('NCBI'),
#     "ITIS": Path('ITIS'),
#     "NCBI_blast": Path('NCBI/blast'),
#     "NCBI_blast_db": Path('NCBI/blast/db'),
#     "NCBI_blast_windowmasker_files": Path('NCBI/blast/windowmasker_files'),
#     "NCBI_pub_taxonomy": Path('NCBI/pub/taxonomy'),
#     "NCBI_refseq_release": Path('NCBI/refseq/release'),
#     "ITIS_taxonomy": Path('ITIS/taxonomy'),
# }

# # Ortholog/utils.py config
# cls = object
# project = "project"
# project_path = "pp"
#
# cls.project = project
# cls.project_path = project_path / Path(project)
# cls.project_index = cls.project_path / Path('index')
# cls.user_index = cls.project_path / Path('index')
#
# cls.user_db = cls.project_path / Path('databases')
# cls.ncbi_db_repo = cls.user_db / Path('NCBI')
# cls.blast_db = cls.ncbi_db_repo / Path('blast') / Path('db')
# cls.windowmaker_files = cls.ncbi_db_repo / Path('blast') / Path('windowmaker_files')
# cls.ncbi_taxonomy = cls.ncbi_db_repo / Path('pub') / Path('taxonomy')
# cls.ncbi_refseq_release = cls.ncbi_db_repo / Path('refseq') / Path('release')
# cls.itis_db_repo = cls.user_db / Path('ITIS')
#
# cls.project_database = cls.user_db / Path(project)
# cls.db_archives = cls.user_db / Path('archive')
#
# cls.raw_data = cls.project_path / Path('raw_data')
# cls.data = cls.project_path / Path('data')
# cls.research_path = cls.project_path

# # config_template_new.yml configuration
#
# """Database_config:
#   email:  "rgilmore@umc.edu"
#   driver: "sqlite3"
#   GenBank_config:
#     download_taxonomy_database:
#       db_type: "biosql"
#       sub_path: "/refseq/release"
#     download_refseq_release_files:
#       collection_subset: "vertebrate_mammalian"
#       seqtype: "rna"
#       filetype: "gbff"
#     upload_refseq_release_files:
#       collection_subset: "vertebrate_mammalian"
#       seqtype: "rna"
#       filetype: "gbff"
#       upload_list: []
#       extension: ".gbk.db"""

# # db_config.yml configuration
#
# """Database_Config:
#   NCBI_blast_db: True
#   NCBI_blast_windowmasker_files: True
#   NCBI_pub_taxonomy: True
#   NCBI_refseq_release: True
#   ITIS_taxonomy: True
#
# Archive_Config:
#   Full: True
#   Projects:
#     flag: False
#     Project_Name_1: False
#     Project_Name_2: False
#     Project_Name_3: False
#   NCBI: False
#   ITIS: False
#   NCBI_blast: False
#   NCBI_blast_db: False
#   NCBI_blast_windowmasker_files: False
#   NCBI_pub_taxonomy: False
#   NCBI_refseq_release: False
#   ITIS_taxonomy: False"""

# Merge No. 3
# Merging of No. 2 and webster_config archive options
# TODO-ROB:  Make changes accordingly in db_config.yml file, cookie_jar.py (bake_the_db_repo; archiving), Ortholog/utils.py (standalone config),
# TODO-ROB: (cont..), config_template_new/existing.yml file (database_management.py conifg_options)

"""
Database_config:
  email:  "rgilmore@umc.edu"
  driver: "sqlite3"
  Full: 
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['']"
  Projects:
    Project_Name_1: 
      configure_flag: True
      archive_flag: True
      delete_flag: False
      path: "!!python/object/apply:pathlib.Path ['Project_Name_1']"
    Project_Name_2: False
      configure_flag: True
      archive_flag: True
      delete_flag: False
      path: "!!python/object/apply:pathlib.Path ['Project_Name_2']"
    Project_Name_3: False
      configure_flag: True
      archive_flag: True
      delete_flag: False
      path: "!!python/object/apply:pathlib.Path ['Project_Name_3']"
  NCBI:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI']"
  ITIS:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['ITIS']"
  NCBI_blast:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI', 'blast']"
  NCBI_blast_db:
    configure_flag: True
    archive_flag:  True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI', 'blast', 'db']"
  NCBI_blast_windowmasker_files:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI', 'blast', 'windowmasker_files']"
  NCBI_pub_taxonomy:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI', 'pub', taxonomy']"
  NCBI_refseq_release:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['NCBI', 'refseq', 'release']"
    upload_flag: True
    db_type: "biosql"
    sub_path: "/refseq/release"
    collection_subset: "vertebrate_mammalian"
    seqtype: "rna"
    filetype: "gbff"
    upload_list: []
    extension: ".gbk.db"
  ITIS_taxonomy:
    configure_flag: True
    archive_flag: True
    delete_flag: False
    path: "!!python/object/apply:pathlib.Path ['ITIS', 'taxonomy']"
"""

# # MERGE No. 2
# # Merging of No. 1 and config_template_new.yml
#
# """Database_config:
#   email:  "rgilmore@umc.edu"
#   driver: "sqlite3"
#   Full:
#     configure_flag: True
#     archive_flag: True
#   Projects:
#     Project_Name_1:
#       configure_flag: True
#       archive_flag: True
#     Project_Name_2: False
#       configure_flag: True
#       archive_flag: True
#     Project_Name_3: False
#       configure_flag: True
#       archive_flag: True
#   NCBI:
#     configure_flag: True
#     archive_flag: True
#   ITIS:
#     configure_flag: True
#     archive_flag: True
#   NCBI_blast:
#     configure_flag: True
#     archive_flag: True
#   NCBI_blast_db:
#     configure_flag: True
#     archive_flag:  True
#   NCBI_blast_windowmasker_files:
#     configure_flag: True
#     archive_flag: True
#   NCBI_pub_taxonomy:
#     configure_flag: True
#     archive_flag: True
#   NCBI_refseq_release:
#     configure_flag: True
#     archive_flag: True
#     upload_flag: True
#     db_type: "biosql"
#     sub_path: "/refseq/release"
#     collection_subset: "vertebrate_mammalian"
#     seqtype: "rna"
#     filetype: "gbff"
#     upload_list: []
#     extension: ".gbk.db"
#   ITIS_taxonomy:
#     configure_flag: True
#     archive_flag: True
# """

# MERGE No. 1
# Merging of "Database_config" and "Archive_config" from db_config.yml file

# """
# Database_Config:
#   Full:
#     configure_flag: True
#     archive_flag: True
#   Projects:
#     Project_Name_1:
#       configure_flag: True
#       archive_flag: True
#     Project_Name_2: False
#       configure_flag: True
#       archive_flag: True
#     Project_Name_3: False
#       configure_flag: True
#       archive_flag: True
#   NCBI:
#     configure_flag: True
#     archive_flag: True
#   ITIS:
#     configure_flag: True
#     archive_flag: True
#   NCBI_blast:
#     configure_flag: True
#     archive_flag: True
#   NCBI_blast_db:
#     configure_flag: True
#     archive_flag:  True
#   NCBI_blast_windowmasker_files:
#     configure_flag: True
#     archive_flag: True
#   NCBI_pub_taxonomy:
#     configure_flag: True
#     archive_flag: True
#   NCBI_refseq_release:
#     configure_flag: True
#     archive_flag: True
#     collection_subset: "vertebrate_mammalian"
#   ITIS_taxonomy:
#     configure_flag: True
#     archive_flag: True
# """