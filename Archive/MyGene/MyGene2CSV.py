import mygene
import pandas as pd

class GetMyGene(object):
    """Use MyGene and a list of refseq_rna accessions to retrieve gene information.

    This script is designed to generate some basic gene information from a list of
    refseqrna accession numbers for human genes.

    The output is a .csv file which includes a gene summary, entrez id, and gene
    full name.

    This pipeline can be used with a number of different input types or 'scopes'
    such as entrezid, gene symbols, etc.
    """
    def __init__(self, filepath):
        """ Create lists for the columns in the file 'Gene', 'Accession'.

        Use mal.keys() to find out the names of the columns.

        Ideally, you should import a csv file with a column titled Accessions with
        a list of refseqrna accessions in order to retrieve the genes.

        More options will be added later.
        """
        self.filepath = filepath
        mal = pd.read_csv(self.filepath)  # List of refseq accessions
        refseq_list = list([accession.upper() for accession in mal.Accessions])
        self.refseq_list = refseq_list

    def retrieve_data(self, outfile, fields='symbol,name,entrezgene,summary', species='human'):
        # Import mygene.MyGeneInfo() search command
        mg = mygene.MyGeneInfo()

        # Create a mygene query handle to get the data
        basic_info = mg.querymany(self.refseq_list, scopes='refseq',
                                  fields=fields,
                                   species=species, returnall=True, as_dataframe=True,
                                   size=1, verbose=True)

        # Use pandas to turn results of the mygene queries into dataframes
        # Reset the index on the dataframe
        basic_info['out'].reset_index(level=0, inplace=True)
        data = basic_info['out']
        gene_info = pd.DataFrame(data)
        gene_info.drop(gene_info.columns[[1,2,6]], axis=1, inplace=True)

        # Rename the columns
        gene_info.rename(columns={'entrezgene': 'Entrez ID','summary':
            'Gene Summary','query': 'RefSeqRNA Accession','name': 'Gene Name'},
            inplace=True)

        # Create the NCBI links using a for loop
        # This is the base url for gene specific pages on NCBI
        baseurl = 'https://www.ncbi.nlm.nih.gov/gene/'

        # Create an empty list that can be appended
        urllist = []

        # Create a for loop that creates the url using the Entrez ID
        # This loop also appends the url to a list and turns it into a link
        for entrezid in gene_info['Entrez ID']:
            entrezid = int(entrezid)
            url = baseurl + str(entrezid)
            # Important step
            # Format the url so that it becomes a hyperlink
            url = '<a href="{0}">{0}</a>'.format(url)
            urllist.append(url)

        # Use pandas to add the urllist & tiers list as columns
        # Turn the ncbiurls list into a dataframe using pandas
        ncbiurls = pd.DataFrame(urllist, columns=['NCBI Link'], dtype=str)

        # List of dataframes I want to combine
        frames = [gene_info, ncbiurls]

        # Merge the dataframes into 1 dataframe
        alldata = pd.concat(frames, axis=1)
#        alldata.rename(columns={'Gene': 'Gene Symbol'}, inplace=True)

#        # Sort the data by the Tier column
#        alldata = alldata.sort_values(['Tier'], ascending=True)

        # Save the merged dataframes to one file that contains all of the information
        alldata.to_csv(outfile, index=False)

