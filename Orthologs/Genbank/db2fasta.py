import os
from BioSQL import BioSeqDatabase
from dir_mana import dir_mana
from lister import Lister

# Set the home directory
home = os.getcwd()

# Set the project name, username, and directory location
project = "GPCR-Orthologs-Project"
user = "rgilmore"
where = dir_mana(home, project)


class Database2Fasta(object):
    """
    Download specific genbank features from a user prefilled BioSQL database
    as fasta files.
    """

    def __init__(self, path='/', fasta_update=False):
        self.FASTA_update = fasta_update
        # Always make sure this file name is correct
        self.what = Lister('MAFV3.1.csv')
        self.org_list = self.what.org_list
        self.gene_list = self.what.gene_list

        if os.path.isdir(where.GB_FASTA + '/Vallender_Data') is False:
            os.mkdir(where.GB_FASTA + '/Vallender_Data')
        self.feature_list = [
            'CDS',
            'Protein',
            'Misc_Features',
            'Other']  # Feature list for directory tree
        path_list = where.path_list_make(where.VERT_MAM, where.NCBI_DATA)
        print('path_list: ', path_list)
        if path == '/':
            self.path, t = self.db2fasta_mkdir(
                where.GB_FASTA + '/Vallender_Data', path_list, self.FASTA_update)
        else:
            self.path = path
            print('self.path: ', self.path)

    def fasta_gather(self):
        """
        Extract the desired fasta feature.
        """
        os.chdir(self.path)
        os.mkdir('MASTER_FASTA')
        os.mkdir('FASTA')
        os.chdir(
            where.VALLENDER_DATA +
            '/refseq/release/multiprocessing/Databases')
        db_name = []
        for file in os.listdir(os.getcwd()):
            if file.endswith('.db'):
                db_name.append(str(file))

        for G_key, G_value in self.what.tier_frame_dict.items():
            Tier = G_key
            db_name = str(Tier) + '.db'
            print('db_name: ', Tier + '.db')
            name = str(db_name)
            server = BioSeqDatabase.open_database(
                driver='sqlite3',
                db=where.VALLENDER_DATA + ('/refseq/release/multiprocessing/Databases/%s' % name))

            os.chdir(self.path + '/FASTA')
            os.mkdir(Tier)
            os.chdir(Tier)
            Tier_path = os.getcwd()
            os.chdir(self.path + '/MASTER_FASTA')
            os.mkdir(Tier)
            os.chdir(Tier)
            M_Tier_path = os.getcwd()
            for Gene in self.what.tier_frame_dict[Tier].T:
                os.chdir(Tier_path)
                os.mkdir(Gene)
                os.chdir(Gene)
                Gene_path = os.getcwd()
                os.chdir(M_Tier_path)
                os.mkdir(Gene)
                os.chdir(Gene)
                M_Gene_path = os.getcwd()
                for Organism in self.org_list:
                    Accession = str(self.what.gene_dict[Gene][Organism])
                    Accession, Sup, Version = Accession.partition('.')
                    Accession = Accession.upper()
                    server_flag = False

                    if server_flag is True:
                        break

                    for sub_db_name in server.keys():

                        try:
                            db = server[sub_db_name]
                            record = db.lookup(accession=Accession)

                            feature_list = []
                            feat_count = 0
                            for feature in record.features:
                                print(feature.type)
                                ref = record.annotations['accessions'][0]
                                GI = record.annotations['gi']
                                feat_count += 1
                                feature_list.append(feature.type)
                                name1 = str(feature.type)

                                # Use the features to initialize
                                if feature.type in feature_list:
                                    x = feature_list.count(feature.type)
                                    name1 = name1 + str(x)
                                if feature.type == "source":
                                    o = feature.qualifiers["organism"]
                                    o[0] = '[' + str(o[0]) + ']'
                                if feature.type == "gene":
                                    c = 0
                                    for x in feature.qualifiers['db_xref']:
                                        c += 1
                                        if 'GI' in x:
                                            AN, sup, GI = x.partition(':')
                                ###############################################
                                """ Change the way the files are saved.  Could just use the copy function.
                                The current way writes a seperate file with the same info"""
                                ###############################################
                                # Create CDS fasta files
                                if feature.type == "CDS":

                                    with open(Gene_path + ('/%s_%s_%s.ffn' % (Gene, Organism, name1)), 'w') as f:
                                        f.write(
                                            ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + str(
                                                feature.extract(record.seq)))
                                        print('File Create: ', f.name)
                                    with open(Gene_path + ('/%s_%s_%s.faa' % (Gene, Organism, 'protein')), 'w') as f:
                                        c = 0
                                        for x in feature.qualifiers['db_xref']:
                                            c += 1
                                            if 'GI' in x:
                                                AN, sup, GI = x.partition(':')
                                        f.write(">" + "gi|" + GI + "|" + "ref" + "|" + str(
                                            feature.qualifiers['protein_id'][0]) + '| ' + str(
                                            feature.qualifiers['product'][0]) + ' ' + o[0] + '\n' + str(
                                            feature.qualifiers['translation'][0]))
                                        print('File Create: ', f.name)

                                    with open(M_Gene_path + ('/MASTER_%s_%s.ffn' % (Gene, name1)), 'a') as f:
                                        f.write(
                                            ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + str(
                                                feature.extract(record.seq)) + "\n\n")
                                        print('File Create: ', f.name)

                                    with open(M_Gene_path + ('/MASTER_%s_%s.faa' % (Gene, 'protein')), 'a') as f:
                                        c = 0
                                        for x in feature.qualifiers['db_xref']:
                                            c += 1
                                            if 'GI' in x:
                                                AN, sup, GI = x.partition(':')
                                        f.write(">" + "gi|" + GI + "|" + "ref" + "|" + str(
                                            feature.qualifiers['protein_id'][0]) + '| ' + str(
                                            feature.qualifiers['product'][0]) + ' ' + o[0] + '\n' + str(
                                            feature.qualifiers['translation'][0]) + '\n\n')
                                        print('File Create: ', f.name)

                                #         ##############################################################################
                                #
                                #
                                #         ##############################################################################
                                #         # Create Misc_Features fasta files
                                # elif feature.type == "misc_feature":
                                #
                                #     with open(home + '/MASTER_Misc_Features/Misc_Features/%s/%s_%s_%s.ffa' % (
                                #     gene_list[line_count - 2], gene_list[line_count - 2], org_list[cell_count - 2], name),
                                #               'w') as f:
                                #         f.write(
                                #             ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + 'Feature: ' + str(
                                #                 feature.qualifiers['note'][0]) + '\n' + str(
                                #                 feature.extract(record.seq)) + "\n")
                                #
                                #     with open(home + '/MASTER_Misc_Features/%s/MASTER_%s_%s.ffa' % (
                                #     gene_list[line_count - 2], gene_list[line_count - 2], n), 'a') as f:
                                #         f.write(
                                #             ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + 'Feature: ' + str(
                                #                 feature.qualifiers['note'][0]) + '\n' + str(
                                #                 feature.extract(record.seq)) + "\n\n")
                                #         ##############################################################################
                                #
                                #
                                #         ##############################################################################
                                #         # Dynamically set up the Other directory tree
                                # elif feature.type != "variation":
                                #     name2 = str(feature.type)
                                #     print(name)
                                #
                                #     if os.path.isdir(home + '/MASTER_Other/MASTER_%s' % name2) is False:
                                #         os.mkdir(home + '/MASTER_Other/MASTER_%s' % name2)
                                #         os.mkdir(home + '/MASTER_Other/MASTER_%s/%s' % (name2, name2))
                                #     if os.path.isdir(home + '/MASTER_Other/MASTER_%s/%s' % (name2, name2)) is False:
                                #         os.mkdir(home + '/MASTER_Other/MASTER_%s/%s' % (name2, name2))
                                #     with open('/work5/r2294/bin/NCBI_data/Input_Files/Gene_List.csv', 'r') as csvfile:
                                #         gene = csv.reader(csvfile, delimiter=',')
                                #         for line in gene:
                                #             for cell in line:
                                #                 if os.path.isdir(home + '/MASTER_Other/MASTER_%s/%s/%s' % (
                                #                 name2, name2, cell)) is False:
                                #                     os.mkdir(home + '/MASTER_Other/MASTER_%s/%s/%s' % (name2, name2, cell))
                                #                 if os.path.isdir(home + '/MASTER_Other/MASTER_%s/%s' % (
                                #                 name2, cell)) is False:
                                #                     os.mkdir(home + '/MASTER_Other/MASTER_%s/%s' % (name2, cell))
                                #                     ##############################################################################
                                #                     # Create Other fasta files
                                #     with open(home + '/MASTER_Other/MASTER_%s/%s/%s/%s_%s_%s.fasta' % (
                                #     name2, name2, gene_list[line_count - 2], gene_list[line_count - 2],
                                #     org_list[cell_count - 2], name1), 'w') as f:
                                #         f.write(
                                #             ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + str(
                                #                 feature.extract(record.seq)) + "\n")
                                #
                                #     with open(home + '/MASTER_Other/MASTER_%s/%s/MASTER_%s_%s.fasta' % (
                                #     name2, gene_list[line_count - 2], gene_list[line_count - 2], n), 'a') as f:
                                #         f.write(
                                #             ">" + "gi|" + GI + "|" + "ref" + "|" + ref + '|' + ' ' + record.description + "\n" + str(
                                # feature.extract(record.seq)) + "\n\n")

                        except IndexError:
                            continue

    def db2fasta_mkdir(self, path, p_list, update):
        if update is True:
            path = where.dir_archive(path, p_list)
        else:
            path = where.dir_make(path, p_list)
        return path
