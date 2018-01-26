"""Helpful utilities for performing Blastn."""
import os
import csv
import time
from datetime import datetime
# import shutil
# import pkg_resources
from importlib import import_module
from multiprocessing.pool import ThreadPool
from pathlib import Path
import pandas as pd
import platform
from warnings import warn
import sqlite3
from sqlalchemy import create_engine

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Tools.otherutils import runcmd

blastutils_log = LogIt().default(logname="blast-utils", logfile=None)
seqidlist_log = LogIt().default(logname="gi-lists", logfile=None)

_datefmt = '%I:%M:%S %p on %m-%d-%Y'
_date = str(datetime.now().strftime(_datefmt))


def map_func(hit):
    """Use the map function for formatting hit id's.
    This will be used later in the script.
    """
    hit.id1 = hit.id.split('|')[3]  # accession number
    hit.id2 = hit.id.split('|')[1]  # gi number
    hit.id = hit.id[:-2]
    return hit


def paml_org_formatter(organisms):
    """Format a list for PAML.

    :param organisms:  Input a list of organisms
    """
    # XXX PAML no longer needs a format different than `Homo_sapiens`
    org_list = []
    for organism in organisms:
        genus, sep, species = organism.partition('_')
        org = ''.join([genus[0], sep, species[0:28]])
        org_list.append(org)
    return org_list


def gene_list_config(file, data_path, gene_list, taxon_dict, logger):
    """Create or use a blast configuration file (accession file).
    This function configures different files for new BLASTS.
    It also helps recognize whether or not a BLAST was terminated in the middle
    of the workflow.  This removes the last line of the accession file if it
    is incomplete.

    :param file:  An accession file to analyze.
    :param data_path:  The path of the accession file.
    :param gene_list:  A gene list in the same order as the accession file.
    :param taxon_dict:  A taxon id dictionary for logging purposes.
    :param logger:  A logger.
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
                org = header[len(row) - 1]  # The last column(organism) accessed in the last row
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
            continued_gene_list = list(x for i, x in enumerate(gene_list, 1) if i > count)
        return continued_gene_list
    # If the file doesn't exist return nothing
    else:
        logger.info("A new BLAST started at %s" % _date)
        return None


def seqid_list_config(seqid_list_path, taxonomy_ids, research_path=None, config=False):
    """Create a seqid list based on the refseq_rna database for each taxonomy id.

    It will also convert the gi list into a binary file which is more
    efficient to use with NCBI's Standalone Blast tools.
    """
    warn("NCBI has deprecated using GI numbers.", OrthoEvolDeprecationWarning)
    if config:
        # Directory and file handling
        raw_data_path = research_path / Path('raw_data')
        index_path = research_path / Path('index')
        taxid_file = index_path / Path('taxids.csv')
        pd.Series(taxonomy_ids).to_csv(str(taxid_file), index=False)

        # TODO Rework this
        create_seqid_lists(seqid_list_path=raw_data_path, taxonomy_ids=taxonomy_ids)

    else:
        create_seqid_lists(seqid_list_path=seqid_list_path, taxonomy_ids=taxonomy_ids)


def create_seqid_lists(seqid_list_path, taxonomy_ids):
    """Use the blastdbcmd tool to generate seqid lists.

    It then uses the blastdb_aliastool to turn the list into a binary file.
    The input (id) for the function is a taxonomy id.
    """
    warn("NCBI has deprecated using GI numbers.", OrthoEvolDeprecationWarning)
    if os.path.exists(str(seqid_list_path)):
        os.chdir(str(seqid_list_path))
        # Use the accession #'s and the blastdbcmd tool to generate gi lists
        # based on Organisms/Taxonomy id's.
        # TODO Create blastdbcmd commandline tools
        gi_time_secs = time.time()
        with ThreadPool(3) as gilist_pool:
            gilist_pool.map(_taxid2seqidlist, taxonomy_ids)
            minutes = (time.time() - gi_time_secs) / 60
        seqidlist_log.info("Took %s minutes to create gi binary files." % minutes)


def _taxid2seqidlist(taxonomy_id):
    """Use a taxonomy id in order to get the list of GI numbers."""
    warn("NCBI has deprecated using GI numbers.", OrthoEvolDeprecationWarning)
    tid = str(taxonomy_id)
    binary = tid + 'gi'

    if binary not in os.listdir():
        if platform.system() == 'Linux':
                # TODO Convert to subprocess
                # TODO Test this on Linux
                runcmd("blastdbcmd -db refseq_rna -entry all -outfmt '%g %a' | awk ' { if ($2 == " + tid + ") { print $1 } } ' > " + tid + "gi.txt")
                seqidlist_log.info(tid + "gi.txt has been created.")

                # Convert the .txt file to a binary file using the blastdb_aliastool
                runcmd("blastdb_aliastool -gi_file_in " + tid + "gi.txt -gi_file_out " + tid + "gi")
                seqidlist_log.info(tid + "gi binary file has been created.")

                # Remove the gi.text file
                os.remove(tid + "gi.txt")
                seqidlist_log.info(tid + "gi.text file has been deleted.")
        else:
            raise NotImplementedError(platform.system() + 'is not supported')
    else:
        seqidlist_log.info('%s already exists' % str(binary))


def my_gene_info(acc_dataframe, blast_query='Homo_sapiens'):
    """Use Biothings' MyGene api to get information about genes.
    :param acc_dataframe:  A pandas dataframe containing the accession csv file data.
    :param blast_query:  The query organism for used during Blasting.
    :return:  Returns a data-frame with hot data about each gene.
    """
    mygene = import_module('mygene')
    blastutils_log.info("Getting Pre-BLAST information about the target genes using MyGene...")
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
                                species='human', returnall=True, as_dataframe=True,
                                size=1, verbose=True)
    # TODO-ROB:  Logging here
    # Turn my-gene queries into a data frame and then reset the index
    mygene_query['out'].reset_index(level=0, inplace=True)
    mg_df = pd.DataFrame(mygene_query['out'])
    mg_df.drop(mg_df.columns[[1, 2, 6]], axis=1, inplace=True)
    # Rename the columns
    mg_df.rename(columns={'entrezgene': 'Entrez ID', 'summary': 'Gene Summary', 'query': 'RefSeqRNA Accession',
                          'name': 'Gene Name'},
                 inplace=True)

    # Create NCBI links using a for loop and the Entrez IDs
    urls = [('<a href="{0}">{0}</a>'.format('https://www.ncbi.nlm.nih.gov/gene/' + str(entrez_id)))
            for entrez_id in mg_df['Entrez ID']]

    ncbi = pd.DataFrame(urls, columns=['NCBI Link'], dtype=str)
    # Merge, sort, and return the my-gene data frame

    hot_data = pd.concat([pd.Series(df.Tier, dtype=str), df.Gene, mg_df, ncbi], axis=1)
    hot_data.rename(columns={'Gene': 'Gene Symbol'}, inplace=True)
    hot_data = hot_data.sort_values(['Tier'], ascending=True)

    return hot_data


def get_dup_acc(acc_dict, gene_list, org_list):
    """
    This function is used to analyze an accession file post-BLAST.  It uses the accession dictionary as a base to get
    duplicated accession numbers.
    :param acc_dict:  A dictionary with accession numbers as keys, and a gene/organism list as values.
    :param gene_list:  A full list of genes.
    :param org_list:  A full list of organisms.
    :return:  A master duplication dictionary used to initialize the duplicate class variables.
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
                    blastutils_log.warning("A duplicate accession number(%s) persists ONLY across %s for %s." % (accession, o, genes))
                    duplicated_dict['organisms'][o][accession] = genes
                    del duplicated_dict['genes'][g]
                    break
                # Duplication across an organisms, but also somewhere else
                elif orgs.count(o) != 1:
                    alt_genes = list(
                        gene for gene, org in go_list if org == o)
                    blastutils_log.warn("A duplicate accession number(%s) persists across %s for %s." % (accession, o, alt_genes))
                    blastutils_log.warn("%s is also duplicated elsewhere." % accession)
                    duplicated_dict['organisms'][o][accession] = alt_genes

                # Duplicates that persist across a gene
                if genes.count(g) == len(go_list):
                    blastutils_log.critical("A duplicate accession number(%s) persists across %s for %s." % (accession, g, orgs))
                    duplicated_dict['genes'][g][accession] = orgs
                    del duplicated_dict['organisms'][o]
                    break
                # Duplication across a gene, but also somewhere else
                elif genes.count(g) != 1:
                    alt_orgs = list(
                        org for gene, org in go_list if gene == g)
                    blastutils_log.critical("A duplicate accession number(%s) persists across %s for %s." % (accession, g, alt_orgs))
                    blastutils_log.critical("%s is also duplicated elsewhere." % accession)
                    duplicated_dict['genes'][g][accession] = alt_orgs

                    # This is the "somewhere else" if the duplication is random or not categorized
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
                    blastutils_log.critical("%s is duplicated, but cannot be categorized as random." % accession)
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


def get_miss_acc(acc_dataframe):
    """Analyze an accession file post BLAST.
    It generates several files and dictionaries regarding missing accession
    numbers.
    :param acc_dataframe: A pandas dataframe containing the accession csv file data(post BLAST).
    :return: A dictionary with data about the missing accession numbers by Gene
             and by Organism.
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
            missing_genes = miss_gene_df.ix[:, organism].to_dict()
            # Do a list comprehension to get a list of genes
            missing_dict['organisms'][organism]['missing genes'] = list(key for key, value in missing_genes.items()
                                                                        if value)  # Value is True for miss accns
            blastutils_log.critical("%s is missing %s." % (organism, str(missing_dict['organisms'][organism]['missing genes'])))
            # Number of missing accessions per organism
            missing_dict['organisms'][organism]['count'] = miss
            total_miss += miss
    missing_dict['organisms']['count'] = total_miss

    # Get missing Accessions by Gene
    gene_dict = miss_gene_df.sum(axis=1).to_dict()
    total_miss = 0
    for gene, miss in gene_dict.items():
        if miss != 0:
            missing_orgs = miss_gene_df.T.ix[:, gene].to_dict()
            missing_dict['genes'][gene] = {}
            # Do a list comprehension to get a list of organisms
            missing_dict['genes'][gene]['missing organisms'] = list(key for key, value in missing_orgs.items()
                                                                    if value  # Value is True for missing accessions
                                                                    if key != 'Tier')  # Don't include 'Tier'
            blastutils_log.critical("%s is missing %s." % (gene, str(missing_dict['genes'][gene]['missing organisms'])))
            # Number of missing accessions per gene
            missing_dict['genes'][gene]['count'] = miss
            total_miss += miss
    missing_dict['genes']['count'] = total_miss

    return missing_dict


def get_pseudogenes():
    """Denote which genes are sudogenes."""
    raise NotImplementedError


def accession_csv2sqlite(acc_file, table_name, db_name, path):
    """
    Convert a OrthoEvolution accession file in csv format
    to an sqlite3 database.

    :param acc_file:  The name of the accession file.  The file name is used to create a table in the
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
    engine = create_engine('sqlite:////%s' % db_path)
    with engine.connect() as conn, conn.begin():
        df = pd.read_csv(acc_path)
        df.to_sql(name=table_name, con=conn, if_exists='append', index=False)


def accession_sqlite2pandas(table_name, db_name, path, exists=True, acc_file=None):
    """
    Convert a sqlite3 database with an OrthoEvolution accession table to a pandas dataframe.
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
        accession_csv2sqlite(acc_file=acc_file, table_name=table_name, db_name=db_name, path=path)
    conn = sqlite3.connect(str(db_path))
    df = pd.read_sql_query("SELECT * FROM %s" % table_name, conn)
    conn.close()
    return df
