# ***********************************************PRE BLAST ANALYSIS TOOLS********************************************* #
# ***********************************************PRE BLAST ANALYSIS TOOLS********************************************* #


def my_gene_info(acc_path, blast_query='Homo_sapiens'):
    import mygene
    import pandas as pd
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
    mg_df.rename(columns={'entrezgene': 'Entrez ID', 'summary':
                         'Gene Summary', 'query': 'RefSeqRNA Accession', 'name': 'Gene Name'},
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
