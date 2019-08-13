    from OrthoEvol.Manager.management import ProjectManagement
    from OrthoEvol.Manager.database_dispatcher import DatabaseManagement
    import yaml

    pm_config_file = "$config_file"

    with open(pm_config_file, 'r') as f:
       pm_config = yaml.load(f)

    pm = ProjectManagement(**pm_config["Management_config"])
    R_R = DatabaseManagement(config_file="$config_file", proj_mana=pm)

    R_R.upload_refseq_release_files(collection_subset="$collection_subset", seqtype="$seqtype", seqformat="$seqformat", upload_list="$sub_list", add_to_default="$add_to_default")