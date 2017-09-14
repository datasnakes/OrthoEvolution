import mygene
import pandas as pd

mg = mygene.MyGeneInfo()


#------------------------------------------------------------------------------
# Create a list of gene symbols/names for .csv file
#------------------------------------------------------------------------------

# Create lists for the 3 columns in the file 'Tier, 'Gene', 'Accession'
# Use tga.keys() to figure out the names of the columns you should use
tg = pd.read_csv('tiers_genes.csv')  # List of refseq
genes_list = list(tg.Gene)
tier_list = list(tg.Tier)

#------------------------------------------------------------------------------
# Set up Input to start command if gene list is correct
#------------------------------------------------------------------------------
"""
x = str(input('Is the input properly formatted? (Type Yes or No) '))
if x == 'Yes':
    print("\n" + "MyGene will start." + "\n")
else:
     raise SystemExit
"""
#------------------------------------------------------------------------------
# Use MyGene to get gene information
#------------------------------------------------------------------------------

# Create a mygene query handle to get the data
basic_gene_info = mg.querymany(genes_list, scopes='symbol',
                          fields='symbol,name,entrezgene,summary',
                           species='human', returnall=True, as_dataframe=True,
                           size=1, verbose=True)


#------------------------------------------------------------------------------
# Use pandas to turn results of the mygene queries into dataframes
#------------------------------------------------------------------------------
"""
Before this step, it'd be wise to check your output to check for duplicates or
missing values.
ex: basic_gene_info.keys()
ex: basic_gene_info['dup']
"""
# Print out the keys and check them
print('The dict keys are ' + str(basic_gene_info.keys()))
print('The output is: ' + str(basic_gene_info['out']))
print('The duplicate results are: ' + str(basic_gene_info['dup']))
print('The missing objects are: ' + str(basic_gene_info['missing']))

# Decide whether to continue with the data as is
input('Continue with using the output? ')

# Turn the dict into a pandas csv file
basic_gene_info['out'].to_csv('basic_gene_info.csv', sep=',', encoding='utf-8')
df = pd.read_csv('basic_gene_info.csv')
data = df
gene_info = pd.DataFrame(data)
gene_info.drop(data.columns[[1,2,6]], axis=1, inplace=True)

# Rename the columns
gene_info.rename(columns={'entrezgene': 'Entrez ID','summary':
    'Gene Summary','query': 'Gene Symbol','name': 'Gene Name'}, inplace=True)

gene_info.to_csv('basic_gene_info.csv', index=False)


#------------------------------------------------------------------------------
# Create the NCBI links using a for loop
#------------------------------------------------------------------------------

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
# Use pandas to add the urllist & tierdata as columns
#------------------------------------------------------------------------------

# Turn the ncbiurls list into a dataframe using pandas
ncbiurls = pd.DataFrame(urllist, columns=['NCBI Link'], dtype=str)
frames = [tier_list, gene_info, ncbiurls]   # List of dataframes I want to combine
alldata = pd.concat(frames, axis=1)  # Merge the dataframes into 1 dataframe
alldata.to_csv('final_gene_info.csv', index=False)