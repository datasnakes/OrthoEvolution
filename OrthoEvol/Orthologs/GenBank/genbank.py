"""Class for managing, downloading and extracting features from genbank files."""
import os
import shutil
from pathlib import Path
from BioSQL import BioSeqDatabase
from Bio import SeqIO

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Orthologs.utils import attribute_config
from OrthoEvol.Orthologs.Blast.orthologs_blastn import OrthoBlastN
from OrthoEvol.Tools.otherutils.other_utils import makedirectory
from OrthoEvol.Orthologs.Blast.comparative_genetics import BaseComparativeGenetics


class GenBank(object):
    """This class will handle GenBank files in various ways."""

    def __init__(self, project, project_path=None, solo=False, multi=True,
                 archive=False, min_fasta=True, blast=OrthoBlastN, **kwargs):
        """Handle GenBank files in various ways.

        It allows for refseq-release .gbff files to be downloaded from NCBI
        and uploaded to a BioSQL database (biopython).  Single .gbk files can be
        downloaded from the .gbff, and uploaded to a custom BopSQL database for
        faster acquisition of GenBank data.

        :param project:  The name of the project.
        :param project_path: The relative path to the project.
        :param solo:  A flag for adding single fasta files.
        :param multi:  A flag for adding multi-fasta files.
        :param archive: A flag for archiving current GenBank Data.  # TODO
        :param min_fasta: A flag for minimizing FASTA file headers.
        :param blast:  The blast parameter is used for composing various
                       Orthologs.Blast classes.  Can be a class, a dict,
                       or none.
        :returns:  .gbff files/databases, .gbk files/databases, & FASTA files.
        """

        # TODO-ROB: Change the way the file systems work.
        self.project = project
        self.solo = solo
        self.multi = multi
        self.min_fasta = min_fasta
        self.genbanklog = LogIt().default(logname="GenBank", logfile=None)

        # Configuration of class attributes
        add_self = attribute_config(self, composer=blast, checker=OrthoBlastN,
                                    checker2=BaseComparativeGenetics,
                                    project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Configuration
        self.target_gbk_db_path = self.user_db / Path(self.project)
        Path.mkdir(self.target_gbk_db_path, parents=True, exist_ok=True)

        # Make a list of BioSQL database(.db) files that contain GenBank info
        self.db_files_list = []
        for FILE in os.listdir(str(self.ncbi_db_repo)):
            if FILE.endswith('.db'):
                self.db_files_list.append(str(FILE))

    @staticmethod
    def name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode):
        """Name a fasta file.

        Provide a uniquely named FASTA file:
        * Coding sequence:
            * Single - "<path>/<gene>_<organism><feat_type_rank>.<extension>"
            * Multi  - "<path>/<gene><feat_type_rank>.<extension>"
        * Other:
            * Single - "<path>/<gene>_<organism>_<feat_type_rank>.<extension>"
            * Multi  - "<path>/<gene>_<feat_type_rank>.<extension>"

        :param path:  The path where the file will be made.
        :param gene:  The gene name.
        :param org:  The organism name.
        :param feat_type:  The type of feature from the GenBank record.
                           (CDS, UTR, misc_feature, variation, etc.)
        :param feat_type_rank:  The feature type  + the rank.
                                (There can be multiple misc_features and
                                 variations)
        :param extension:  The file extension.
                           (".ffn", ".faa", ".fna", ".fasta")
        :param mode:  The mode ("w" or "a") for writing the file.  Write to a
                      solo-FASTA file.  Append a multi-FASTA file.
        :return:  The uniquely named FASTA file.
        """

        # Create path variables.  (typically raw_data/<gene>/GENBANK
        feat_path = path
        # Create a format-able string for file names
        if feat_type_rank is "CDS":
            single = '%s_%s%s%s'
            multi = '%s%s%s'
        else:
            single = '%s_%s_%s%s'
            multi = '%s_%s%s'
        # Create different names based on fasta file type
        if mode == 'w':
            file_path = feat_path / Path(single % (gene, org, feat_type_rank, extension))
        elif mode == 'a':
            file_path = feat_path / Path(multi % (gene, feat_type_rank, extension))

        # Make the base directory and return an open file.
        makedirectory(feat_path)
        file = open(file_path, mode)
        return file

    @staticmethod
    def protein_gi_fetch(feature):
        """Retrieve the protein gi number.

        :param feature:  Search the protein feature for the GI number.
        :return:  The protein GI number as a string.
        """

        # Find the protein gi number under the features qualifiers.
        for x in feature.qualifiers:
            if 'GI' in x:
                _, _, p_gi = x.partition(':')
                return p_gi

    def create_post_blast_gbk_records(self, org_list, gene_dict):
        """Create a single GenBank file for each ortholog.

        After a blast has completed and the accession numbers have been compiled
        into an accession file, this class searches a local NCBI refseq release
        database composed of GenBank records.  This method will create a single
        GenBank file (.gbk) for each ortholog with an accession number.
        The create_post_blast_gbk_records is only callable if the
        the instance is composed by one of the Blast classes.  This method also
        requires an NCBI refseq release database to be set up with the proper
        GenBank Flat Files (.gbff) files.

        :param org_list:  List of organisms
        :param gene_dict:  A nested dictionary for accessing accession numbers.
                           (e.g. gene_dict[GENE][ORGANISM} yields an accession
                           number)
        :return:  Does not return an object, but creates genbank files.
        """

        # Parse the tier_frame_dict to get the tier
        for G_KEY, _ in self.tier_frame_dict.items():
            tier = G_KEY
            # Parse the tier based transformed dataframe to get the gene
            for GENE in self.tier_frame_dict[tier].T:
                # Parse the organism list to get the desired accession number
                for ORGANISM in org_list:
                    accession = str(gene_dict[GENE][ORGANISM])
                    accession, _, version = accession.partition('.')
                    # When parsing a GenBank database, the version needs to be removed.
                    accession = accession.upper()
                    server_flag = False
                    # Search the databases and create a GenBank file.
                    self.get_gbk_file(accession, GENE, ORGANISM, server_flag=server_flag)

    def get_gbk_file(self, accession, gene, organism, server_flag=None):
        """Search a GenBank database for a target accession number.

        This function searches through the given NCBI databases (created by
        uploading NCBI refseq .gbff files to a BioPython BioSQL database) and
        creates single GenBank files.  This function can be used after a
        blast or on its own.  If used on it's own then the NCBI .db files must
        be manually moved to the proper directories.

        :param accession: Accession number of interest without the version.
        :param gene: Target gene of the accession number parameter.
        :param organism: Target organism of the accession number parameter.
        :param server_flag:  (Default value = None)
        :return:
        """

        gene_path = self.raw_data / Path(gene) / Path('GENBANK')
        Path.mkdir(gene_path, parents=True, exist_ok=True)

        # Parse each database to find the proper GenBank record
        for FILE in self.db_files_list:
            db_file_path = self.ncbi_db_repo / Path(FILE)
            # Stop searching if the GenBank record has been created.
            if server_flag is True:
                break
            server = BioSeqDatabase.open_database(driver='sqlite3',
                                                  db=str(db_file_path))
            # Parse the sub-databases
            for SUB_DB_NAME in server.keys():
                db = server[SUB_DB_NAME]
                try:
                    record = db.lookup(accession=accession)
                    gbk_file = '%s_%s.gbk' % (gene , organism)
                    gbk_file_path = gene_path / Path(gbk_file)
                    with open(gbk_file_path, 'w') as GB_file:
                        GB_file.write(record.format('genbank'))
                        self.genbanklog.info(GB_file.name, 'created')
                    # Make sure we have the correct GenBank file.
                    self.gbk_quality_control(gbk_file_path, gene, organism)
                    # Stop searching if the GenBank record has been created.
                    server_flag = True
                    break
                except IndexError:
                    self.genbanklog.critical('Index Error in %s.  Moving to the next database...' % SUB_DB_NAME)
                    continue

        # If the file has not been created after searching, then raise an error
        if server_flag is not True:
            self.genbanklog.critical("The GenBank file was not created for %s (%s, %s)." % (accession, gene, organism))
            raise FileNotFoundError

    def gbk_quality_control(self, gbk_file, gene, organism):
        """Ensures the quality or validity of the retrieved genbank record.

    It takes the GenBank record and check to make sure the Gene and Organism
        from the GenBank record match the Gene and Organism from the accession
        file.  If not, then the Blast has returned the wrong accession number.

        :param gbk_file:  The path to a GenBank file.
        :param gene:  A gene name from the Accession file.
        :param organism:  A gene name from the Accession file.
        :return:
        """

        # TODO-ROB:  Check the bad data here against the misssing/duplicate files
        record = SeqIO.read(gbk_file, 'genbank')
        gene_flag = False
        organism_flag = False
        accession = record.id
        self.gbk_gene_synonym = {}
        self.duplicated_dict["validated"] = {}
        # Get the organism name from the GenBank file
        gbk_organism = record.features[0].qualifiers["organism"]  # A list with one entry
        if len(gbk_organism) == 1:
            gbk_organism = gbk_organism[0]
            gbk_organism = gbk_organism.replace(" ", "_")
        else:
            self.genbanklog.critical("Two organisms exist in the GenBank file.  Is this normal?")
            raise BrokenPipeError

        # Check to make sure the organism in the GenBank file matches the
        # organism from the accession file
        if gbk_organism == organism:
            self.genbanklog.info("The GenBank organism, %s, has been verified for %s." % (organism, gene))
        else:
            organism_flag = True

        # Get the gene from the GenBank files
        gbk_genes = record.features[1].qualifiers["gene"]
        # Get the synonyms from the GenBank file if they exist and add them to
        # the list.
        if "gene_synonym" in str(record.features[1].qualifiers.keys()):
            base_gene_name = gbk_genes
            gbk_genes.extend(record.features[1].qualifiers["gene_synonym"])
            # Create a dictionary from the synonyms
            self.gbk_gene_synonym[base_gene_name] = []
            self.gbk_gene_synonym[base_gene_name].extend(gbk_genes)

        # Check to make sure the gene in the GenBank file matches the gene from
        # the accession file
        for gbk_gene in gbk_genes:
            if gbk_gene == gene:
                gene_flag = False
                self.genbanklog.info("The GenBank gene, %s, has been verified for %s." % (gene, organism))
                break
            else:
                gene_flag = True

        # TODO-ROB:  Add a verified key to the duplicates dictionary.
        # Raise errors.
        if organism_flag is True and gene_flag is True:
            self.genbanklog.critical("The organisms don't match.\n\tGenBank: %s \n\tAccession File: %s" %
                                     (gbk_organism, organism))
            self.genbanklog.critical("The genes don't match. \n\tGenBank: %s \n\tAccession File: %s" %
                                     (gbk_genes, gene))
            raise BrokenPipeError
        elif organism_flag is True:
            self.genbanklog.critical("The organisms don't match.\n\tGenBank: %s \n\tAccession File: %s" %
                                     (gbk_organism, organism))
            raise BrokenPipeError
        elif gene_flag is True:
            self.genbanklog.critical("The genes don't match. \n\tGenBank: %s \n\tAccession File: %s" %
                                     (gbk_genes, gene))
            raise BrokenPipeError

        self.duplicated_dict["validated"][accession] = [gene, organism]

    def gbk_upload(self):
        """Upload a BioSQL database with target GenBank data (.gbk files).

        This method is only usable after creating GenBank records with this
        class.  It uploads a BioSQL databases with target GenBank data (.gbk
        files).  This creates a compact set of data for each project.

        :return:  Does not return an object.
        """

        t_count = 0
        # Parse the tier dictionary
        for TIER in self.tier_frame_dict.keys():
            db_name = str(TIER) + '.db'
            db_file_path = self.target_gbk_db_path / Path(db_name)
            # Create the db file if it exists
            if os.path.isfile(str(db_file_path)) is False:
                self.genbanklog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))

            # If it already exists then the database is bad, or needs to be update.  Delete it.
            else:
                # TODO-ROB: This part is broken until the template db creation and management is added
                os.remove(str(db_file_path))
                self.genbanklog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))

            server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
            gene_path = self.raw_data
            # Parse the raw_data folder to get the name of each gene.
            for GENE in os.listdir(str(gene_path)):
                sub_db_name = GENE
                genbank_path = gene_path / Path(GENE) / Path('GENBANK')
                # Parse the GenBank file names for each gene in order to upload them to a custom BioSQL database
                for FILE in os.listdir(str(genbank_path)):
                    # Try to load the database.
                    try:
                        if sub_db_name not in server.keys():
                            server.new_database(sub_db_name)
                        db = server[sub_db_name]
                        count = db.load(SeqIO.parse(FILE, 'genbank'))
                        server.commit()
                        self.genbanklog.info('Server Commited %s' % sub_db_name)
                        self.genbanklog.info('%s database loaded with %s.' % (db.dbid, FILE))
                        self.genbanklog.info("That file contains %s genbank records." % str(count))
                        t_count = t_count + count
                        self.genbanklog.info('The total number of files loaded so far is %i.' % t_count)
                    # If the database cannot be loaded then rollback the server and raise an error.
                    except BaseException:
                        server.rollback()
                        # Try to delete the sub database and commit
                        try:
                            del server[sub_db_name]
                            server.commit()
                        # If it cannot be deleted then raise an error.
                        except BaseException:
                            raise
                        raise

    def get_fasta_files(self, acc_dict, db=True):
        """Create FASTA files for each GenBank record in the accession dictionary.

        It can search through a BioSQL database or it can crawl a directory
        for .gbk files.

        :param acc_dict:  An accession dictionary like the one created by
                          CompGenObjects.
        :param db:  A flag that determines whether or not to use the custom
                    BioSQL database or to use .gbk files.
                    (Default value = True)
        :return:  Returns FASTA files for each GenBank record.
        """

        # Get FASTA files from the BioSQL GenBank databases.
        if db is True:
            # Parse the directory that contains the databases for the project of interest.
            for database in os.listdir(str(self.target_gbk_db_path)):
                server = BioSeqDatabase.open_database(driver="sqlite3", db=database)
                try:
                    for db_name in server.keys():
                        db = server[db_name]
                        # For each GenBank record in the database write a set of FASTA files.
                        for item in db.keys():
                            record = db.lookup(item)
                            self.write_fasta_files(record, acc_dict)
                            self.genbanklog.info("FASTA files for %s created from BioSQL database." % item)
                except:
                    raise()
        # Get FASTA files from the GenBank files.
        # TODO-ROB change this.  Broken by new directory structure
        # TODO-ROB directory looks like /raw_data/Gene_1/GENBANK/*.gbk
        elif db is False:
            # Parse the directory that contain the GenBank records for the project of interest.
            for _, _, gbk_files in os.walk(str(self.target_gbk_files_path)):
                # For each genbank record write a set of FASTA files.
                for gbk_file in gbk_files:
                    if Path(gbk_file).suffix is '.gbk':
                        record = SeqIO.read(gbk_file, 'genbank')
                        self.write_fasta_files(record, acc_dict)
                        self.genbanklog.info("FASTA files for %s created." % gbk_file)

    def write_fasta_files(self, record, acc_dict):
        """Create a dictionary for formatting the FASTA header & sequence.

        :param record:  A GenBank record created by BioPython.
        :param acc_dict:  Accession dictionary from the CompGenObjects class.
        :return:
        """

        feat_type_list = []
        for feature in record.features:
            # XXX Set up variables to use for dictionary values !!!
            # Basic variables.
            accession = record.id
            gene = acc_dict[accession][0]
            organism = acc_dict[accession][1]
            # Variable for minimalistic FASTA files.
            genus, sep, species = organism.partition('_')
            min_org = str(''.join([genus[0], sep, species[0:28]]))

            # Keep a list of feature types to identify duplicates (for naming the FASTA files).
            # The first iteration of the feature type contains no number.
            # The following iterations are concatenated with a number.
            feat_type = str(feature.type)
            feat_type_list.append(feat_type)
            duplicate_num = feat_type_list.count(feat_type)
            if duplicate_num == 1:
                feat_type_rank = feat_type
            else:
                feat_type_rank = feat_type + str(duplicate_num)
            # XXX END !!!

            # TODO-ROB:  Remove the GI number stuff here or at least prepare for
            # file with no GI.
            # Create a dictionary and format FASTA file entries.
            fmt = {
                'na_gi': str(record.annotations['gi']),
                'aa_gi': str(self.protein_gi_fetch(feature)),
                'na_acc_n': str(accession),
                'aa_acc_n': str(feature.qualifiers['protein_id'][0]),
                'na_description': str(record.description),
                'aa_description': str(feature.qualifiers['product'][0]),
                'na_seq': str(feature.extract(record.seq)),
                'aa_seq': str(feature.qualifiers['translation'][0]),
                'na_misc_feat': str(feature.qualifiers['note'][0]),
                'org': str(organism),
                'gene': str(gene),
                'min_org': str(min_org),
                'feat_type': str(feat_type),
                'feat_type_rank': str(feat_type_rank),
                'path': str(self.raw_data / Path(gene) / Path('GENBANK'))
            }
            # Set up minimalistic FASTA headers and sequence entries for Nucleic Acid and Amino Acid sequences.
            na_entry = ">{min_org}\n{na_seq}\n".format(**fmt)
            aa_entry = ">{min_org}\n{aa_seq}\n".format(**fmt)
            # For full FASTA headers/sequences set min_fasta to False
            if self.min_fasta is False:
                na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description}\n{na_seq}\n".format(**fmt)
                aa_entry = ">gi|{aa_gi}|reg|{aa_acc_n}| {aa_description} {org}\n{aa_seq}\n".format(**fmt)
            # ######### End ######### #

            # ############ Write desired FASTA files ############ #
            if self.solo is True:
                self.solo_fasta(na_entry, aa_entry, fmt)
            if self.multi is True:
                self.multi_fasta(na_entry, aa_entry, fmt)

    def solo_fasta(self, na_entry, aa_entry, fmt):
        """This method writes a sequence of a feature to a uniquely named file using a dictionary for formatting.

        :param na_entry:  A string representing the Nucleic Acid sequence data in FASTA format.
        :param aa_entry:  A string representing the Amino Acid sequence data in FASTA format.
        :param fmt:  A dictionary for formatting the FASTA entries and the file names.
        :return:  Does not return an object, but creates single entry FASTA files.
        """

        mode = 'w'

        # Create the desired variables from the formatter dictionary.
        feat_type = fmt['feat_type']
        feat_type_rank = fmt['feat_type_rank']
        path = fmt['path']
        gene = fmt['gene']
        org = fmt['org']

        if feat_type == "CDS":
            # Create a .ffn file (FASTA for Coding Nucleic Acids)
            extension = '.ffn'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
            # Create a .faa file (FASTA for Amino Acids)
            extension = '.faa'
            file = self.name_fasta_file(path, gene, org, 'Protein', feat_type_rank, extension, mode)
            file.write(aa_entry)
            file.close()

        elif feat_type == "misc_feature":
            # Create a custom entry for miscellaneous features.
            na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description} Feature: {na_misc_feat}\n{na_seq}\n".format(**fmt)
            # Creates .fna files (generic FASTA file for Nucleic Acids)
            extension = '.fna'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()

        elif feat_type != "variation":
            # Creates .fasta files (generic FASTA file)
            extension = '.fasta'
            file = self.name_fasta_file(path, gene, org, 'Other', feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()

    def multi_fasta(self, na_entry, aa_entry, fmt):
        """Append an othologous sequence of a feature to a uniquely named file.

        Usese a dictionary for formatting.

        :param na_entry:  A string representing the Nucleic Acid sequence data in FASTA format.
        :param aa_entry:  A string representing the Amino Acid sequence data in FASTA format.
        :param fmt:  A dictionary for formatting the FASTA entries and the file names.
        :return:  Does not return an object, but creates or appends to a multi entry FASTA file.
        """

        mode = 'a'

        # Create the desired variables from the formatter dictionary.
        feat_type = fmt['feat_type']
        feat_type_rank = fmt['feat_type_rank']
        path = fmt['path']
        gene = fmt['gene']
        org = fmt['org']

        if feat_type == "CDS":
            # Create a MASTER .ffn file (multi-FASTA file for Coding Nucleic Acids)
            extension = '.ffn'
            file = self.name_fasta_file(path, gene, org, feat_type,
                                        feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
            # Create a MASTER .faa file (multi-FASTA file for Amino Acids)
            extension = '.faa'
            file = self.name_fasta_file(path, gene, org, feat_type,
                                        feat_type_rank, extension, mode)
            file.write(aa_entry)
            file.close()
        elif feat_type == "misc_feature":
            na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description} Feature: {na_misc_feat}\n{na_seq}\n".format(**fmt)
            # Creates .fna files (generic FASTA file for Nucleic Acids)
            extension = '.fna'
            file = self.name_fasta_file(path, gene, org, feat_type,
                                        feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
