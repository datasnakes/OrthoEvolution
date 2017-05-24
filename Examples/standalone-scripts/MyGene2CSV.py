"""
"""
# Modules Used
import mygene
import pandas as pd

#------------------------------------------------------------------------------
# Create lists for the 3 columns in the file 'Tier, 'Gene', 'Homo_Sapiens'
# Use tga.keys() to figure out the names of the columns you should use
mal = pd.read_csv('Master_Accession_File_Final.csv')  # List of refseq
refseq_list = list([accession.upper() for accession in mal.Homo_sapiens])

# Must be uppercase for mygene to recognize
#refseq_list = [accession.upper() for accession in refseq_list]

#------------------------------------------------------------------------------
# Import mygene.MyGeneInfo() search command
mg = mygene.MyGeneInfo()

# Create a mygene query handle to get the data
basic_info = mg.querymany(refseq_list, scopes='refseq',
                          fields='symbol,name,entrezgene,summary',
                          species='human', returnall=True, as_dataframe=True,
                          size=1, verbose=True)

#------------------------------------------------------------------------------
# Use pandas to turn results of the mygene queries into dataframes
# Reset the index on the dataframe
basic_info['out'].reset_index(level=0, inplace=True)
data = basic_info['out']
gene_info = pd.DataFrame(data)
gene_info.drop(gene_info.columns[[1, 2, 6]], axis=1, inplace=True)

# Rename the columns
gene_info.rename(columns={'entrezgene': 'Entrez ID', 'summary':
                          'Gene Summary', 'query': 'RefSeqRNA Accession', 'name': 'Gene Name'},
                 inplace=True)

#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Use pandas to add the urllist & tiers list as columns
# Turn the ncbiurls list into a dataframe using pandas
ncbiurls = pd.DataFrame(urllist, columns=['NCBI Link'], dtype=str)

# List of dataframes I want to combine
frames = [mal.Tier, mal.Gene, gene_info, ncbiurls]

# Merge the dataframes into 1 dataframe
alldata = pd.concat(frames, axis=1)
alldata.rename(columns={'Gene': 'Gene Symbol'}, inplace=True)

# Sort the data by the Tier column
alldata = alldata.sort_values(['Tier'], ascending=True)

# Save the merged dataframes to one file that contains all of the information
alldata.to_csv('datasets/gpcr_gene_info.csv', index=False)
