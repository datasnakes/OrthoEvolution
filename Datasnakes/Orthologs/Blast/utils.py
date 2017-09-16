"""Helpful utilities for performing Blastn"""
import pandas as pd
import os
import csv
from pathlib import Path
import time
from importlib import import_module


def map_func(hit):
    """Use the map function for formatting hit id's.
    This will be used later in the script.
    """
    hit.id1 = hit.id.split('|')[3]  # accession number
    hit.id2 = hit.id.split('|')[1]  # gi number
    hit.id = hit.id[:-2]
    return hit


# XXX PAML no longer needs a format different than `Homo_sapiens`
def paml_org_formatter(organisms):
    org_list = []
    for organism in organisms:
        genus, sep, species = organism.partition('_')
        org = ''.join([genus[0], sep, species[0:28]])
        org_list.append(org)
    return org_list


def gene_list_config(file, data_path, gene_list, taxon_dict, logger):
    """Create or use a blast configuration file.
    This function configures different files for new BLASTS.
    It also helps recognize whether or not a BLAST was terminated
    in the middle of the dataset.  This removes the last line of
    the accession file if it is incomplete.
    """

    ending = gene = org = taxid = None
    output_dir_list = os.listdir(data_path)  # Make a list of files
    # If the file exists then make a gene list that picks up from the last BLAST
    if file in output_dir_list:
        header = pd.read_csv(str(Path(data_path) / Path(file)), dtype=str)
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

            # ######### Start logging ######### #
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
            # ######## End Logging ######## #
            # The continued gene list starts with the previous gene.
            continued_gene_list = list(x for i, x in enumerate(gene_list, 1) if i > count)
        return continued_gene_list
    # If the file doesn't exist return nothing
    else:
        logger.info("A new BLAST started at %s" % time.time())
        return None

def my_gene_info(acc_path, blast_query='Homo_sapiens'):
    mygene = import_module('mygene')

    # Initialize variables and import my-gene search command
    urls = []
    df = pd.read_csv(str(acc_path), dtype=str)
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
    """Get duplicate accessions.
    This function is used to analyze an accession file post-BLAST.
    It uses the accession dictionary as a base.

    :return: A master duplication dictionary used to initialize the
    duplicate class variables.
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
                    duplicated_dict['organisms'][o][accession] = genes
                    del duplicated_dict['genes'][g]
                    break
                # Duplication across an organisms, but also somewhere else
                elif orgs.count(o) != 1:
                    alt_genes = list(
                        gene for gene, org in go_list if org == o)
                    duplicated_dict['organisms'][o][accession] = alt_genes

                # Duplicates that persist across a gene
                if genes.count(g) == len(go_list):
                    duplicated_dict['genes'][g][accession] = orgs
                    del duplicated_dict['organisms'][o]
                    break
                # Duplication across a gene, but also somewhere else
                elif genes.count(g) != 1:
                    alt_orgs = list(
                        org for gene, org in go_list if gene == g)
                    duplicated_dict['genes'][g][accession] = alt_orgs

                    # This is the "somewhere else" if the duplication is random or not categorized
                    # The duplication is random
                if genes.count(g) == 1 and orgs.count(o) == 1:
                    del duplicated_dict['organisms'][o]
                    del duplicated_dict['genes'][g]
                    if accession not in duplicated_dict['random']:
                        duplicated_dict['random'][accession] = []
                    duplicated_dict['random'][accession].append(go)
                    # There is another category of duplication that I'm missing
                    # TODO-ROB:  If an other exists throw a warning in the
                    # logs
                else:
                    del duplicated_dict['organisms'][o]
                    del duplicated_dict['genes'][g]
                    if accession not in duplicated_dict['other']:
                        duplicated_dict['other'][accession] = []
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


def get_miss_acc(acc_file_path):
    """
    This function is used to analyze an accession file post BLAST.
    It generates several files and dictionaries regarding missing accession numbers.

    :param acc_file_path: An accession file (post BLAST).
    :return: A dictionary with data about the missing accession numbers by Gene and
    by Organism.
    """
    # TODO-ROB: Add Entrez ID validation;  Get info from xml files???
    missing_dict = dict()
    missing_dict['organisms'] = {}
    missing_dict['genes'] = {}

    # Initialize the Data Frames
    data = pd.read_csv(str(acc_file_path), dtype=str)
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

            # Number of missing accessions per gene
            missing_dict['genes'][gene]['count'] = miss
            total_miss += miss
    missing_dict['genes']['count'] = total_miss

    return missing_dict


def get_pseudogenes():
    """ UNDER DEVELOPMENT!!!
    This subclass will denote which genes are sudogenes.
    """
    print(__doc__)