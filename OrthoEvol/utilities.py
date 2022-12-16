"""Helpful utilities for performing Blastn."""
# Standard Library
import csv
import itertools
import os
import shutil
import sqlite3
import subprocess as sp
import sys
import pkg_resources
from threading import Timer
from datetime import datetime
from importlib import import_module
from pathlib import Path
from subprocess import TimeoutExpired
from tempfile import TemporaryFile
# BioPython
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
# OrthoEvol
from OrthoEvol.Cookies.cookie_jar import Oven
from OrthoEvol.Tools.logit import LogIt
# Other
import pandas as pd
import yaml

# Set up logging
blastutils_log = LogIt().default(logname="blast-utils", logfile=None)
seqidlist_log = LogIt().default(logname="gi-lists", logfile=None)
utils_log = LogIt().default(logname="utils", logfile=None)
_datefmt = '%I:%M:%S %p on %m-%d-%Y'
_date = str(datetime.now().strftime(_datefmt))


class BlastUtils(object):

    def __init__(self):
        """
        Various utilities to help with blast specific functionality.
        """
        pass

    def map_func(self, hit):
        """Format/parse hit ids generated from blast xml results.
        """
        hit.id1 = hit.id.split('|')[3]  # accession number
        hit.id2 = hit.id.split('|')[1]  # gi number
        hit.id = hit.id[:-2]
        return hit

    def paml_org_formatter(self, organisms):
        """Take a list of organisms and format each organism name for PAML,
        which can only take names that are less than a certain length
        (36 characters?).

        :param organisms:  A list of organisms
        :type organisms:  list.
        """
        # XXX PAML no longer needs a format different than `Homo_sapiens`
        org_list = []
        for organism in organisms:
            genus, sep, species = organism.partition('_')
            org = ''.join([genus[0], sep, species[0:28]])
            org_list.append(org)
        return org_list

    def gene_list_config(self, file, data_path, gene_list, taxon_dict, logger):
        """Create or use a blast configuration file (accession file).

        This function configures different files for new BLASTS.
        It also helps recognize whether or not a BLAST was terminated in the
        middle of the workflow.  This removes the last line of the accession
        file if it is incomplete.

        :param file:  An accession file to analyze.
        :type file:  str.
        :param data_path:  The path of the accession file.
        :type data_path:  str.
        :param gene_list:  A gene list in the same order as the accession file.
        :type gene_list:  list.
        :param taxon_dict:  A taxon id dictionary for logging purposes.
        :type taxon_dict:  dict.
        :param logger:  A LogIt logger for logging.
        :type logger:  LogIt.
        :return:  Returns a continued gene_list to pick up from an interrupted Blast.
        """

        ending = gene = org = taxid = None
        output_dir_list = os.listdir(str(data_path))  # Make a list of files

        # If the file exists, make a gene list that picks up from the last BLAST
        if file in output_dir_list:
            header = pd.read_csv(os.path.join(data_path, file), dtype=str)
            header = header.axes[1].tolist()
            file = Path(data_path) / Path(file)
            with open(str(file), 'r') as open_file:
                csv_file = csv.reader(open_file)
                count = 0
                # Iterate the csv file and determine the last gene to be blasted
                for row in csv_file:
                    if count == 0:
                        header = row
                    count += 1
                    ending = row  # The last row
                    gene = ending[1]  # The last row's gene
                    # The last column(organism) accessed in the last row
                    org = header[len(row) - 1]
                    taxid = taxon_dict[org]  # The taxon id of the organism

                # Start logging
                ncbi = str("""result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna", strand="plus",
                    evalue=0.001, out="%s_%s.xml", outfmt=5, gilist=%s + "gi", max_target_seqs=10, task="blastn")"""
                           % (gene, org, taxid))
                logger.warning("An incomplete accession file was produced from the previous BLAST, which was terminated "
                               "midway through the procedure.")
                logger.info("The last row looks like: \n\t%s\n\t%s\n" % (header, ending))
                logger.info("The BLAST ended on the following query: \n%s" % ncbi)
                if len(ending) < len(header):
                    logger.info("Restarting the BLAST for the previous gene...")
                    count = count - 2
                # End logging
                # The continued gene list starts with the previous gene.
                continued_gene_list = list(
                    x for i, x in enumerate(
                        gene_list, 1) if i > count)
            return continued_gene_list
        # If the file doesn't exist return nothing
        else:
            logger.info("A new BLAST started at %s" % _date)
            return None

    def my_gene_info(self, acc_dataframe, blast_query='Homo_sapiens'):
        """Use Biothings' MyGene api to get information about genes.

        :param acc_dataframe:  A pandas dataframe containing the accession csv file data.
        :type acc_dataframe:  pd.DataFrame.
        :param blast_query:  The query organism for used during Blasting.
        :type blast_query:  str.
        :return:  Returns a data-frame with hot data about each gene.
        :rtype:  pd.DataFrame.
        """

        mygene = import_module('mygene')
        blastutils_log.info(
            "Getting Pre-BLAST information about the target genes using MyGene...")
        # Initialize variables and import my-gene search command
        urls = []
        df = acc_dataframe
        df_gene = df.set_index('Gene')
        blast_query_list = df_gene[blast_query].tolist()
        mg = mygene.MyGeneInfo()

        # Create a my-gene query handle to get the data
        human = list(x.upper() for x in blast_query_list)
        mygene_query = mg.querymany(human, scopes='refseq',
                                    fields='symbol,name,entrezgene,summary',
                                    species='human', returnall=True,
                                    as_dataframe=True,
                                    size=1, verbose=True)
        # TODO-ROB:  Logging here
        # Turn my-gene queries into a data frame and then reset the index
        mygene_query['out'].reset_index(level=0, inplace=True)
        mg_df = pd.DataFrame(mygene_query['out'])
        mg_df.drop(mg_df.columns[[1, 2, 6]], axis=1, inplace=True)
        # Rename the columns
        mg_df.rename(columns={'entrezgene': 'Entrez ID',
                              'summary': 'Gene Summary',
                              'query': 'RefSeqRNA Accession',
                              'name': 'Gene Name'},
                     inplace=True)

        # Create NCBI links using a for loop and the Entrez IDs
        urls = [('<a href="{0}">{0}</a>'.format('https://www.ncbi.nlm.nih.gov/gene/' + str(entrez_id)))
                for entrez_id in mg_df['Entrez ID']]

        ncbi = pd.DataFrame(urls, columns=['NCBI Link'], dtype=str)
        # Merge, sort, and return the my-gene data frame

        query_data = pd.concat([pd.Series(df.Tier, dtype=str),
                                df.Gene, mg_df, ncbi], axis=1)
        query_data.rename(columns={'Gene': 'Gene Symbol'}, inplace=True)
        query_data = query_data.sort_values(['Tier'], ascending=True)

        return query_data

    def get_dup_acc(self, acc_dict, gene_list, org_list):
        """Get duplicated accession numbers during post-blast analysis.

        :param acc_dict:  A dictionary with accession numbers as keys, and a
                          gene/organism list as values.
        :type acc_dict:  dict.
        :param gene_list:  A full list of genes.
        :type gene_list:  list.
        :param org_list:  A full list of organisms.
        :type org_list:  list.
        :return:  A master duplication dictionary used to initialize the
                  duplicate class variables.
        :rtype:  dict.
        """

        duplicated_dict = dict()
        duplicated_dict['accessions'] = {}
        duplicated_dict['genes'] = {}
        duplicated_dict['organisms'] = {}
        duplicated_dict['random'] = {}
        duplicated_dict['other'] = {}
        acc_dict = acc_dict
        dup_gene_count = {}
        dup_org_count = {}
        dup_acc_count = {}

        for accession, go_list in acc_dict.items():
            # Finding duplicates by using the length of the accession
            # dictionary
            if len(go_list) > 1:
                # dict['accessions']['XM_000000'] = [[g, o], [g, o]]
                duplicated_dict['accessions'][accession] = go_list
                genes, orgs = zip(*go_list)
                genes = list(genes)
                orgs = list(orgs)
                # Process the duplicates by categorizing and storing in a
                # dictionary
                for go in go_list:
                    g = go[0]
                    o = go[1]
                    # Initialize the dictionaries if they haven't already been
                    if g not in duplicated_dict['genes']:
                        duplicated_dict['genes'][g] = {}
                    if o not in duplicated_dict['organisms']:
                        duplicated_dict['organisms'][o] = {}
                        # Categorize the different types of duplication
                    # Duplicates that persist across an organisms
                    if orgs.count(o) == len(go_list):
                        blastutils_log.warning(
                            "A duplicate accession number(%s) persists ONLY across %s for %s." % (accession, o, genes))
                        duplicated_dict['organisms'][o][accession] = genes
                        del duplicated_dict['genes'][g]
                        break
                    # Duplication across an organisms, but also somewhere else
                    elif orgs.count(o) != 1:
                        alt_genes = list(
                            gene for gene, org in go_list if org == o)
                        blastutils_log.warning(
                            "A duplicate accession number(%s) persists across %s for %s." % (accession, o, alt_genes))
                        blastutils_log.warning(
                            "%s is also duplicated elsewhere." % accession)
                        duplicated_dict['organisms'][o][accession] = alt_genes

                    # Duplicates that persist across a gene
                    if genes.count(g) == len(go_list):
                        blastutils_log.critical(
                            "A duplicate accession number(%s) persists across %s for %s." % (accession, g, orgs))
                        duplicated_dict['genes'][g][accession] = orgs
                        del duplicated_dict['organisms'][o]
                        break
                    # Duplication across a gene, but also somewhere else
                    elif genes.count(g) != 1:
                        alt_orgs = list(
                            org for gene, org in go_list if gene == g)
                        blastutils_log.critical(
                            "A duplicate accession number(%s) persists across %s for %s." % (accession, g, alt_orgs))
                        blastutils_log.critical(
                            "%s is also duplicated elsewhere." % accession)
                        duplicated_dict['genes'][g][accession] = alt_orgs

                        # This is the "somewhere else" if the duplication
                        # is random or not categorized
                        # The duplication is random
                    if genes.count(g) == 1 and orgs.count(o) == 1:
                        del duplicated_dict['organisms'][o]
                        del duplicated_dict['genes'][g]
                        if accession not in duplicated_dict['random']:
                            duplicated_dict['random'][accession] = []
                        blastutils_log.critical("%s is randomly duplicated." % accession)
                        duplicated_dict['random'][accession].append(go)
                        # There is another category of duplication that I'm missing
                        # TODO-ROB:  If an other exists throw a warning in the
                        # logs
                    else:
                        del duplicated_dict['organisms'][o]
                        del duplicated_dict['genes'][g]
                        if accession not in duplicated_dict['other']:
                            duplicated_dict['other'][accession] = []
                        blastutils_log.critical(
                            "%s is duplicated, but cannot be categorized as random." % accession)
                        duplicated_dict['other'][accession].append(go)
            # Duplicate Organism count dictionary
            dup_org = pd.DataFrame.from_dict(duplicated_dict['organisms'])
            for org in org_list:
                try:
                    dup_org_count[org] = dup_org[org].count()
                except KeyError:
                    dup_org_count[org] = 0
            # Duplicate Gene count dictionary
            dup_gene = pd.DataFrame.from_dict(duplicated_dict['genes'])
            for gene in gene_list:
                try:
                    dup_gene_count[gene] = dup_gene[gene].count()
                except KeyError:
                    dup_gene_count[gene] = 0
            # Duplicate Accession count dictionary
            for accn, go in duplicated_dict['accessions'].items():
                dup_acc_count[accn] = go.__len__()
        return duplicated_dict

    def get_miss_acc(self, acc_dataframe):
        """Get missing accession numbers during post-blast analysis.

        :param acc_dataframe:  A pandas dataframe containing the accession csv file data(post BLAST).
        :type acc_dataframe:  pd.DataFrame
        :return:  A dictionary with data about the missing accession numbers by Gene and by Organism.
        :rtype:  dict.
        """
        # TODO-ROB: Add Entrez ID validation;  Get info from xml files???
        missing_dict = dict()
        missing_dict['organisms'] = {}
        missing_dict['genes'] = {}

        # Initialize the Data Frames
        data = acc_dataframe
        df = data.set_index('Gene')
        miss_gene_df = df.isnull()
        miss_org_df = df.T.isnull()

        # Get missing Accessions by Organism
        organism_dict = miss_org_df.sum(axis=1).to_dict()
        total_miss = 0
        for organism, miss in organism_dict.items():
            if miss != 0:
                missing_dict['organisms'][organism] = {}
                # Missing Gene dict {'HTR1A': True}
                missing_genes = miss_gene_df.loc[:, organism].to_dict()
                # Do a list comprehension to get a list of genes
                missing_dict['organisms'][organism]['missing genes'] = list(key for key, value in missing_genes.items()
                                                                            if value)  # Value is True for miss accns
                blastutils_log.critical("%s is missing %s." % (organism, str(
                    missing_dict['organisms'][organism]['missing genes'])))
                # Number of missing accessions per organism
                missing_dict['organisms'][organism]['count'] = miss
                total_miss += miss
        missing_dict['organisms']['count'] = total_miss

        # Get missing Accessions by Gene
        gene_dict = miss_gene_df.sum(axis=1).to_dict()
        total_miss = 0
        for gene, miss in gene_dict.items():
            if miss != 0:
                missing_orgs = miss_gene_df.T.loc[:, gene].to_dict()
                missing_dict['genes'][gene] = {}
                # Do a list comprehension to get a list of organisms
                missing_dict['genes'][gene]['missing organisms'] = list(key for key, value in missing_orgs.items()
                                                                        if value  # Value is True for missing accessions
                                                                        if key != 'Tier')  # Don't include 'Tier'
                blastutils_log.critical("%s is missing %s." % (
                    gene, str(missing_dict['genes'][gene]['missing organisms'])))
                # Number of missing accessions per gene
                missing_dict['genes'][gene]['count'] = miss
                total_miss += miss
        missing_dict['genes']['count'] = total_miss

        return missing_dict

    # def get_pseudogenes(self):
    #     """Denote which genes are pseudogenes."""
    #     raise NotImplementedError

    def accession_csv2sqlite(self, acc_file, table_name, db_name, path):
        """Convert am OrthoEvol csv accession file to an sqlite3 database.

        :param acc_file:  The name of the accession file.  The file name is
                          used to create a table in the
        sqlite3 database.  Any periods will be replaced with underscores.
        :type acc_file: str
        :param table_name: The name of the table in the database.
        :type table_name: str
        :param db_name: The name of the new database.
        :type db_name: str
        :param path: The relative path of the csv file and the database.
        :type path: str
        """
        acc_path = Path(path) / Path(acc_file)
        db_path = Path(path) / Path(db_name)
        with sqlite3.connect(str(db_path)) as conn:
            df = pd.read_csv(acc_path, dtype=str)
            df.to_sql(name=table_name, con=conn, if_exists='replace',
                      index=False)

    def accession_sqlite2pandas(self, table_name, db_name, path, exists=True,
                                acc_file=None):
        """Convert a sqlite3 database with an OrthoEvol accession table to a pandas dataframe.

        :param table_name: Name of the table in the database.
        :type table_name: str
        :param db_name: The name of the new database.
        :type db_name: str
        :param path: The relative path of the csv file and the database.
        :type path: str
        :param exists: A flag used to create a database if needed.
        :type exists: bool
        :param acc_file:  The name of the accession file.  The file name is used to create a table in the
        sqlite3 database.  Any periods will be replaced with underscores.
        :type acc_file: str
        :return:
        :rtype:
        """
        db_path = Path(path) / Path(db_name)
        if not exists or not db_path.is_file():
            self.accession_csv2sqlite(acc_file=acc_file, table_name=table_name,
                                      db_name=db_name, path=path)
        conn = sqlite3.connect(str(db_path))
        df = pd.read_sql_query("SELECT * FROM %s" % table_name, conn)
        conn.close()
        return df


class GenbankUtils(object):

    def __init__(self):
        """Various utilities to help with genbank specific functionality."""
        pass

    def multi_fasta_manipulator(self, target_file, man_file, output_file,
                                manipulation='remove'):
        # Inspired by the BioPython Tutorial and Cookbook ("20.1.1 Filtering a
        # sequence file")
        """Manipulate reference sequences in multifasta files.

        The original purpose was to filter files created by the GUIDANCE2
        alignment program, but this function has been made in order to
        accommodate other situations as well.

        :param target_file:  Target multi-FASTA file.
        :type target_file:  str.
        :param man_file:  Sequences in a multi-FASTA file used for manipulation.
        :type man_file:  str.
        :param manipulation:  The type of manipulation.  (remove, add, sort, tbd..)
        :type manipulation:  str.
        :param output_file:  The name of the output multi-FASTA file.
        :type output_file:  str.
        :return:  A multi-FASTA file that has been manipulated accordingly.
        :rtype:  str.
        """
        # Create path variables
        new_file = Path(target_file).parent / Path(output_file)
        # Create a new multi-fasta record object using the target_file,
        # reference, and output
        # Remove specific sequences from a fasta file
        if manipulation == 'remove':
            self.multi_fasta_remove(target_file, man_file, new_file)
        # Combine all the FASTA sequence in one record object
        elif manipulation == 'add':
            self.muli_fasta_add(target_file, man_file, new_file)
        # Sort one fasta file based on the order of another
        # Works for alignments
        elif manipulation == 'sort':
            self.multi_fasta_sort(target_file, man_file, new_file)

        utils_log.info('A new fasta file has been created.')
        return new_file

    # def dir_config(path, tier_frame_dict):
    #     """
    #     Configure the genbank directories.
    #     :param path: Path to create directory structure.
    #     :param tier_frame_dict:  Dictionary from the blast super class.
    #     :return:  Creates a directory structure as follows
    #         --Tier_1
    #             --Gene_1
    #             --Gene_M
    #         --Tier_N
    #             --Gene_M+1
    #             --Gene_N
    #     """
    #     for G_KEY in tier_frame_dict.keys():
    #         tier = G_KEY
    #         tier_path = path / Path(tier)
    #         Path.mkdir(tier_path, parents=True, exist_ok=True)
    #         for GENE in tier_frame_dict[tier].T:
    #             gene_path = tier_path / Path(GENE)
    #             Path.mkdir(gene_path)
    def multi_fasta_remove(self, target_file, man_file, output_file):
        """Remove selected reference sequences in a multi-FASTA files.

        :param target_file:  Target multi-FASTA file.
        :type target_file:  str.
        :param man_file:  A multi-FASTA file with sequences used for removal from the the target file.
        :type man_file:  str.
        :param output_file:  The name of the output multi-FASTA file.
        :type output_file:  str.
        :return:  A multi-FASTA file with removed sequences.
        :rtype:  str.
        """
        rem_file = output_file.stem + '_removed' + output_file.suffix
        rem_file = output_file.parent / Path(rem_file)
        # Turn the reference_file into set of ids
        if os.path.isfile(man_file):
            ids = set(record.id for record in SeqIO.parse(man_file, 'fasta'))
        elif isinstance(man_file, list):
            ids = man_file

        new_records = (record for record in SeqIO.parse(
            target_file, 'fasta') if record.id not in ids)
        old_records = (
            record for record in SeqIO.parse(
                target_file,
                'fasta') if record.id in ids)

        utils_log.info('Sequences have been filtered.')
        SeqIO.write(new_records, str(output_file), 'fasta')
        SeqIO.write(old_records, str(rem_file), 'fasta')

    def muli_fasta_add(self, target_file, man_file, output_file):
        """Add selected reference sequences in a multi-FASTA files.

        :param target_file:  Target multi-FASTA file.
        :type target_file:  str.
        :param man_file:  A multi-FASTA file with sequences used for appending to the target file.
        :type man_file:  str.
        :param output_file:  The name of the output multi-FASTA file.
        :type output_file:  str.
        :return:  A multi-FASTA file with added sequences.
        :rtype:  str.
        """
        # TODO-ROB:  Check for duplicates.
        # Concatenate the multifasta files together by chaining the SeqIO.parse generators
        # Allows one to overwrite a file by using temporary files for storage
        # adding generators -
        # https://stackoverflow.com/questions/3211041/how-to-join-two-generators-in-python
        if os.path.isfile(man_file):
            with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
                new_records = itertools.chain(SeqIO.parse(
                    target_file, 'fasta', ), SeqIO.parse(man_file, 'fasta'))
                count = SeqIO.write(new_records, tmp_file, 'fasta')
                tmp_file.seek(0)
                utils_log.info('temp file count: ' + str(count))
                SeqIO.write(SeqIO.parse(tmp_file, 'fasta'), str(output_file), 'fasta')
            utils_log.info('Sequences have been added.')
        else:
            utils_log.warning('You can only add files together.  Not python objects.')

    def multi_fasta_sort(self, target_file, man_file, output_file):
        """Sorts selected reference sequences in a multi-FASTA files.

        :param target_file:  Target multi-FASTA file.
        :type target_file:  str.
        :param man_file:  A multi-FASTA file with sequences used for sorting the target file.
        :type man_file:  str.
        :param output_file:  The name of the output multi-FASTA file.
        :type output_file:  str.
        :return:  A multi-FASTA file with sorted sequences.
        :rtype:  str.
        """
        # TODO-ROB:  Check for duplicates.
        with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
            aln = MultipleSeqAlignment([])

        # For a reference fasta file make a tuple from ids
        if os.path.isfile(man_file):
            sorted_list = []
            for s in SeqIO.parse(man_file, 'fasta'):
                sorted_list.append(s.id)
            sorted_handle = tuple(sorted_list)
        # For a reference list port in as a tuple
        elif isinstance(man_file, list):
            sorted_handle = tuple(man_file)

        # Parse the reference tuple above the unsorted file
        for sorted_seq_record in sorted_handle:
            for unsorted_aln_record in SeqIO.parse(target_file, 'fasta'):
                # If an appropriate id is found, then append to the MSA object.
                if unsorted_aln_record.id == sorted_seq_record.id:
                    utils_log.info(unsorted_aln_record.id)
                    utils_log.info(sorted_seq_record.id)
                    aln.append(unsorted_aln_record)  # MSA object
                    break
        count = AlignIO.write(aln, tmp_file, 'fasta')
        tmp_file.seek(0)
        utils_log.info('temp file count: ' + str(count))
        AlignIO.write(AlignIO.read(tmp_file, 'fasta'), str(output_file), 'fasta')
        utils_log.info('Alignment has been sorted.')


class OrthologUtils(BlastUtils, GenbankUtils):

    def __init__(self):
        """Various utilities to help with ortholog specific functionality."""
        BlastUtils.__init__(self)
        GenbankUtils.__init__(self)

    def attribute_config(self, cls, composer, checker, project=None,
                         project_path=None, checker2=None):
        """Set or configure attributes.

        Attribute Configuration takes an instance of a class and sets various
        attributes. The attributes are set by determining the type of configuration.
        The class attributes can be composed by another class, they can be set with
        a dictionary, or they can be set using the basic project template.

        :param cls: An instance of a class that will retain the attributes.
        :type cls:  object.
        :param composer: A class that will yield attributes to the cls parameter.
        :type composer:  object.
        :param checker: A checker class used to check the type of the composer.
                Dictionary composers will be treated differently.
        :type checker:  object.
        :type checker:  dict.
        :param project:  The name of the project.
        :type project:  str.
        :param project_path:  The relative path of the project.
        :type project_path:  str.
        :return:  Returns the instance (cls) with new attributes.
        :rtype:  object.
        """
        ac_log = LogIt().default(logname="%s" % cls.__class__.__name__, logfile=None)
        if checker2:
            check2 = issubclass(type(composer), checker2)
        else:
            check2 = None
        # Attribute configuration using checker composition.
        if issubclass(type(composer), checker) or check2:
            for key, value in composer.__dict__.items():
                setattr(cls, key, value)
            ac_log.info("The attribute configuration was accomplished by composing %s with %s." % (
                cls.__class__.__name__, composer.__class__.__name__))

        # Attribute configuration using a dictionary.
        elif isinstance(composer, dict):
            for key, value in composer.items():
                setattr(cls, key, value)
            ac_log.info(
                "The attribute configuration of %s was accomplished by using a dictionary." % cls.__class__.__name__)

        # Attribute configuration without composer
        elif composer is None:
            if not (project or project_path):
                raise BrokenPipeError(
                    "Without the Project Management class, a project name and project path must be included.")
            cls = self.standalone_config(cls, project, project_path)
            ac_log.info(
                "The attribute configuration of %s was accomplished by using a standalone project." % cls.__class__.__name__)
        # Make sure self.project and self.project_path have values
        if not (cls.project or cls.project_path):
            raise BrokenPipeError(
                "The project name and project path attributes have not been set.")

        return cls

    def standalone_config(self, cls, project, project_path, custom=None):
        """Configure a standalone project.

        A standalone configuration uses the variables listed.  These variables are
        either mapped to a basic project, used in a custom configuration, or they
        are mapped to some basic project directories with some custom options.

        :param cls: An instance of a class that will retain the attributes.
        :type cls:  object.
        :param project: The name of the project.
        :type project:  str.
        :param project_path: The relative path of a project.
        :type project_path:  str.
        :param new: The new project flag.
        :type new:  bool.
        :param custom: The custom flag which can be None or a dictionary.
        :type custom:  dict.
        :return: Returns the instance (cls) with new attributes.
        :rtype:  object.
        """

        if not project or not project_path:
            raise NameError('A project and project_path must be specified.')
        cls.project = project
        cls.project_path = Path(project_path) / Path(project)
        cls.project_index = cls.project_path / Path('index')
        cls.user_index = cls.project_path / Path('index')
        cls.db_archives = cls.project_path / Path('archive')
        cls.raw_data = cls.project_path / Path('raw_data')
        cls.data = cls.project_path / Path('data')
        cls.research_path = cls.project_path
        cls.user_db = cls.project_path / Path('databases')
        cls.project_database = cls.user_db / Path(cls.project)
        cls.itis_db_repo = cls.user_db / Path('ITIS')
        cls.ncbi_db_repo = cls.user_db / Path('NCBI')
        cls.blast_db = cls.ncbi_db_repo / Path('blast') / Path('db')
        cls.windowmaker_files = cls.ncbi_db_repo / \
            Path('blast') / Path('windowmaker_files')
        cls.ncbi_taxonomy = cls.ncbi_db_repo / Path('pub') / Path('taxonomy')
        cls.NCBI_refseq_release = cls.ncbi_db_repo / Path('refseq') / Path('release')

        # Use the basic_project cookie to create the directory structure
        Kitchen = Oven(project=project, basic_project=True)
        Kitchen.bake_the_project(cookie_jar=project_path)

        # Use the custom dictionary to set the path variables
        # and to make the directories if necessary.  This overrides
        if custom:
            for key, value in custom.items():
                setattr(cls, key, value)
                if not Path(str(value)).is_dir():
                    Path.mkdir(value, exist_ok=True)

        return cls


class ManagerUtils(object):
    def __init__(self):
        """Various utilities to help with management specific functionality."""
        pass

    def parse_db_config_file(self, config_file):
        """Parse a YAML config file and return the config strategies and keyword arguments.

        :param config_file:  A YAML config file for database management.
        :type config_file:   str.
        :return:  The database config strategies, and the key word arguments for database management.
        :rtype:  tuble.
        """
        kw = {}
        db_config_strategy = {}
        with open(config_file, 'r') as cf:
            db_config = yaml.load(cf, Loader=yaml.FullLoader)
            # Get the configuration for the desired strategy
            for key, value in db_config["Database_config"].items():
                if isinstance(value, dict):
                    db_config_strategy[key] = value
                    continue
                # Get the parameters for the Base class
                else:
                    kw[key] = value
        return db_config_strategy, kw

    def refseq_jobber(self, email_address, base_jobname, id, code, activate, config_dict):
        """Submit python code as a string.

        :param email_address:  The email address for PBS job notification.
        :type email_address:  str.
        :param base_jobname:  The base job name used for the PBS job.  Contains a %s for string formatting.
        :type base_jobname:  str.
        :param id:  An id used to format the base_jobname.
        :type id:  int.
        :param code:  Python code as a string.
        :type code:  str.
        :param activate:  The path to the activate script for the virtual environment being used in the PBS job.
        :type activate:  str.
        """
        job = SGEJob(email_address=email_address, base_jobname=base_jobname % str(id), activate=activate,
                     config=config_dict)
        job.submit_pycode(code=code, wait=False, cleanup=False)

    # def template_jobber(email_address, base_jobname, id, code):
    #     job = SGEJob(email_address=email_address, base_jobname=base_jobname)
    #     job.submit_pycode(code=code, wait=True, cleanup=True)


class CookieUtils(object):
    def __init__(self):
        """Various utilities to help with cookie specific functionality."""
        self.archive_options = {
            "Full": Path(''),
            "NCBI": Path('NCBI'),
            "ITIS": Path('ITIS'),
            "NCBI_blast": Path('NCBI/blast'),
            "NCBI_blast_db": Path('NCBI/blast/db'),
            "NCBI_blast_windowmasker_files": Path('NCBI/blast/windowmasker_files'),
            "NCBI_pub_taxonomy": Path('NCBI/pub/taxonomy'),
            "NCBI_refseq_release": Path('NCBI/refseq/release'),
            "ITIS_taxonomy": Path('ITIS/taxonomy'),
        }

        self.bytesize_options = {
            "B": 1,
            "KB": 1024,
            "MB": 1048576,
            "GB": 1073741824,
            "TB": 1099511627776
        }

    def archive(self, database_path, archive_path, option, delete_flag=False):
        """Archive a database directory from a Cookie templated directory structure.

        This utility creates a YAML config dictionary that contains path-like
        objects for archiving.  The original data
        can be moved to the archive path or deleted all together.

        :param database_path:  A path to a folder that consists of the desired data.
        :type database_path:  str.
        :param archive_path:  A path to an output folder for archived data.
        :type archive_path:  str.
        :param option:  An option for the archiving strategy.  Will be one of the keys in the archive_options.
        :type option:  str.
        :param delete_flag:  A flag for deleting the original data.  USE WITH CAUTION.
        :type delete_flag:  bool.
        :return:  Returns a list of paths to the *.tar.xz archive of the data and/or a path to the original data.
        :rtype:  list.
        """
        archive_dict = {}
        archive_list = []
        archive_log = LogIt().default(logname="Archive", logfile=None)

        if option == "Full":
            full_path = Path(database_path) / self.archive_options["Full"]
            for folder in os.listdir(str(full_path)):
                if os.path.isdir(folder):
                    archive_dict[folder] = database_path / Path(folder)
        elif isinstance(option, list):
            for opt in option:
                other_path = Path(database_path) / self.archive_options[opt]
                archive_dict[opt] = other_path
        else:
            other_path = Path(database_path) / self.archive_options[option]
            archive_dict[option] = other_path

        for arch_name, data_path in archive_dict.items():
            root_dir = str(data_path.parent)
            base_dir = str(data_path.stem)
            d = datetime.now().strftime(format="%Y-%m-%d_%H%M")
            output_pathname = archive_path / Path(arch_name + "." + d)
            # Archive the desired data.
            data_size = self.get_size(start_path=str(data_path))
            archive_log.info("Archiving %s of data." % data_size)
            archive_filename = shutil.make_archive(base_name=str(output_pathname), format="xztar", root_dir=root_dir,
                                                   base_dir=base_dir, logger=archive_log)
            archive_size = self.get_size(archive_filename)
            archive_log.warning("A %s archive file was created at %s." %
                                (archive_filename, archive_size))
            # TODO-ROB:  Logging.  And log to a README.md file.
            # Delete the files if desired.
            if delete_flag:
                archive_log.critical(
                    "The original data will be deleted recursively at %s." % data_path)
                from OrthoEvol import OrthoEvolWarning
                OrthoEvolWarning(
                    "You're about to delete your database (%s).  Are you sure??" % data_path)
                shutil.rmtree(path=data_path)
                archive_list.append(str(archive_filename))
            else:
                archive_log.critical("The original data will be moved recursively from %s to %s." % (
                    data_path, output_pathname))
                output_pathname.mkdir()
                shutil.move(src=str(data_path), dst=str(output_pathname))
                shutil.move(src=str(archive_filename), dst=str(output_pathname))
                archive_list.append(str(output_pathname))

            Path(data_path).mkdir(parents=True, exist_ok=True)
        return archive_list

    def get_size(self, start_path, units="KB"):
        """Determine the size of a directory or a file with the desired units.

        :param start_path:  A file or path for sizing up.
        :type start_path:  str.
        :param units:  The denomination of bytes to return.
        :type units:  str.
        :return:  The size as a string.  (e.g. "3.6 KB")
        :rtype:  str.
        """
        total_size = 0
        if os.path.isfile(start_path):
            size = os.path.getsize(start_path)
            size = str(size / self.bytesize_options[units]) + (" %s" % units)
            return size

        for dirpath, _, filenames in os.walk(start_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total_size += os.path.getsize(fp)
        total_size = str(total_size / self.bytesize_options[units]) + (" %s" % units)
        return total_size


class PackageVersion(object):
    """Get the version of an installed python package."""

    def __init__(self, packagename):
        self.packagename = packagename
        self._getversion()

    def _getversion(self):
        import_module(self.packagename)
        version = pkg_resources.get_distribution(self.packagename).version
        utils_log.info('Version %s of %s is installed.' % (version, self.packagename))


class FunctionRepeater(object):
    """Repeats a function every interval. Ref: https://tinyurl.com/yckgv8m2"""

    def __init__(self, interval, function, *args, **kwargs):
        self._timer = None
        self.function = function
        self.interval = interval
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False


class FullUtilities(CookieUtils, ManagerUtils, OrthologUtils):

    def __init__(self):
        """A class to help bring utility functions to other modules."""
        CookieUtils.__init__(self)
        ManagerUtils.__init__(self)
        OrthologUtils.__init__(self)

    def system_cmd(self, cmd, timeout=None, print_flag=True,
                   write_flag=False, file_name=None, **kwargs):
        """
        A function for making system calls, while preforming proper exception handling.

        :param cmd:  A list that contains the arguments for Popen.
        :param timeout:  A timeout variable.
        :param write_flag:  A flag used to write each line of output to a file.
        :param file_name:  The name of the to be written to.
        :param kwargs:  Any other keyword arguments for Popen.
        :return:  Returns the stdout and stderr as a tuple.
        """
        proc = sp.Popen(cmd, **kwargs, encoding="utf-8")
        stdout_list = []
        for line in iter(proc.stdout.readline, ""):
            stdout_list.append(line.rstrip())
            if print_flag:
                utils_log.info(line, end="")
            sys.stdout.flush()
            if write_flag:
                with open(file_name, 'a') as output_file:
                    output_file.write(line)
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
            if proc.returncode != 0:
                raise sp.CalledProcessError(proc.returncode, cmd, stdout, stderr)
        except TimeoutExpired:
            proc.kill()
        proc.stdout = stdout_list
        return proc

    def group_files_by_size(self, file_dict, _path=None, groups=8):
        """Create a list of dictionaries that contain the specified number of groups of files.

        Each group of files has a total size that are close to each other.

        :param file_dict:  The file names are keys, and file sizes are the values.
        :type file_dict:  dict.
        :param _path:  If the file_dict isn't given, then a path can be given for generating the file_dict.
        :type _path:  str.
        :param groups:  The number of groups to split the files into.
        :type groups:  int.
        :return:  A list with {groups} number of dictionaries.
        :rtype:  list.
        """
        # Generate a file dict from the path if file_dict isn't given.
        if not file_dict:
            file_dict1 = {}
            for file in os.listdir(_path):
                if '.db' in file:
                    continue
                file_dict1[file] = Path(file).stat().st_size
        # Initialize variables
        group_max_size = sum(list(file_dict.values())) / groups
        avg_file_size = sum(list(file_dict.values())) / len(file_dict.keys())
        total_files = len(list(file_dict.values()))
        group_list = []
        file_count = 0
        # Begin parsing the dictionary
        for _ in range(0, groups):
            group_dict = {}
            group_size = 0
            for file in dict(file_dict).keys():
                # Don't add files to the group if it's over the max size.
                if group_size < group_max_size:
                    file_count = file_count + 1
                    group_size = group_size + file_dict[file]
                    group_diff = group_size - group_max_size
                    if (group_size < group_max_size) or (group_diff < avg_file_size):
                        group_dict[file] = file_dict.pop(file)
                # Append to the group list.
                else:
                    group_list.append(group_dict)
                    group_dict = {}
                    break
            # Append the last group to the list.  It gets skipped above.
            if len(group_dict.keys()) > 0:
                group_list.append(group_dict)
        # If there are extra files for some reason add them to one of the groups.
        if file_count != total_files:
            group_list[1].update(file_dict)
        return group_list
