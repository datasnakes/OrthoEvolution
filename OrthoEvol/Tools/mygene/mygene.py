"""Get gene information from the mygene api at http://mygene.info/ """
import mygene
import pandas as pd

from OrthoEvol.Tools.logit import LogIt
# TODO add ability for custom fields
# TODO add ability for custom species.


class MyGene(object):
    """Import a csv of refseq accessions & get gene information from mygene."""

    def __init__(self, infile, outfile):
        """Initialize my gene handle and refseq/accessions list.

        Get the basic gene information. It's best to use a csv file and title
        the row of the accessions list `Accessions`.

        :param infile: Input a csv file with a `Homo_sapiens` column
        :param outfile: Desired path & name of the output file.
        """
        self.infile = infile
        self.outfile = outfile

        self.mygene_log = LogIt().default('mygene', None)

        self.mg = mygene.MyGeneInfo()  # Set up mygene handle
        self.accessions_list = self._import_accfile()  # Create accessions list

        self.fields = 'symbol,name,entrezgene,summary'  # Default fields
        self.species = 'human'  # Species to use.

    def _import_accfile(self):
        """Import the accession file and make Homo_Sapiens column a list.

        :return: Returns a list of accessions.
        """

        accfile = pd.read_csv(self.infile)
        acclist = list([accession.upper() for accession in accfile.Homo_Sapiens])
        return acclist

    def query_mygene(self):
        """Query mygene for gene information."""

        self.mygene_log.info('You are querying: %s' % self.accessions_list)
        basic_info = self.mg.querymany(self.accessions_list, scopes='refseq',
                                       fields=self.fields,
                                       species=self.species, returnall=True,
                                       as_dataframe=True, size=1, verbose=True)

        # basic_info['out'] is the output dataframe.
        # Use basic_info.keys() to find dict keys
        # Reset the index on the dataframe so that each column is on the same
        # level
        basic_info['out'].reset_index(level=0, inplace=True)
        data = basic_info['out']
        gene_info = pd.DataFrame(data)
        gene_info.drop(gene_info.columns[[0, 1, 2]], axis=1, inplace=True)
        gene_info.rename(columns={'symbol': 'Gene Symbol',
                                  'entrezgene': 'Entrez ID',
                                  'name': 'Gene Name',
                                  'summary': 'Summary'}, inplace=True)

        # Create the NCBI links using a for loop
        baseurl = 'https://www.ncbi.nlm.nih.gov/gene/'

        # Create an empty list that can be appended to
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

        # Turn the ncbiurls list into a dataframe using pandas
        ncbiurls = pd.DataFrame(urllist, columns=['NCBI Link'], dtype=str)

        # List of dataframes I want to combine
        frames = [gene_info, ncbiurls]

        # Merge the dataframes into 1 dataframe
        alldata = pd.concat(frames, axis=1)

        # Save the merged dataframes to a file
        alldata.to_csv(self.outfile, index=False)
        self.mygene_log.info('%s has been created.' % str(self.outfile))
