# -*- coding: utf-8 -*-
"""
Integrate MyGene with the Gene web

"""
#------------------------------------------------------------------------------
# Modules Used
#------------------------------------------------------------------------------
import mygene
import csv
import pandas as pd
import sys
from bokeh.charts import HeatMap, show, output_file
from bokeh.palettes import RdYlGn4
from bokeh.models import HoverTool, Legend

# Import myge.MyGeneInfo() search command
mg = mygene.MyGeneInfo()


#------------------------------------------------------------------------------
# Create a list of gene symbols/names for .csv file
#------------------------------------------------------------------------------
g = open('genes.csv')  # 1st column - gene names
genes_list = []   # Initialize a list of genes
genes_list.append('')
file2 = csv.reader(g)
for gene in file2:    # Format a list of genes
    genes = str(gene)
    genes = genes.replace("'", "")
    genes = genes.replace("[", "")
    genes = genes.replace("]", "")
    genes = genes.replace(" ", "_")
    genes_list.append(genes)
print(genes_list)

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
"""
Call querymany method.
Input is "symbol", and you want "entrezgene" (Entrez gene ids) back.

Set as_dataframe to True will return a pandas dataframe object
Set verbose to False as this will suppress the messages like "finished".
The resuls will be a list of dictionaries.
The dictionary contains the entrezid for the "entrezgene" field.
If you want the ensembl ids, use fields='ensembl.gene'

List of fields: http://mygene.info/metadata/fields

Examples:
entrez_ids = mg.querymany(genes_list, scopes='symbol', fields='entrezgene',
                          species='human', returnall=True, as_dataframe=True)

ensembl_ids = mg.querymany(genes_list, scopes='symbol', fields='ensembl.gene',
                          species='human', returnall=True, as_dataframe=True)
"""
# This creates a dictionary of basic human gene information to be used later
basic_info = mg.querymany(genes_list, scopes='symbol',
                          fields='symbol,name,entrezgene,summary',
                           species='human', returnall=True, as_dataframe=True,
                           size=1)

#------------------------------------------------------------------------------
# Use pandas to turn results of the mygene queries into dataframes
#------------------------------------------------------------------------------
"""
Write the dataframe to stdout.
It will not have NaNs on the screen if no matches were found.
Save the data as a .csv file.
Use df.drop to delete columns of the data you don't want.

Additional dictionary command:
To return a dictionary of MyGene.info metadata, use metadata = mg.metadata
"""
# Turn the dict into a pandas csv file
basic_info['out'].to_csv('basic_gene_info.csv', sep=',', encoding='utf-8')
df = pd.read_csv('basic_gene_info.csv')
data = df
gene_info = pd.DataFrame(data)
gene_info.drop(data.columns[[1,2,6]], axis=1, inplace=True)

# Rename the columns
gene_info.rename(columns={'entrezgene': 'Entrez ID','summary':
    'Gene Summary','query': 'Gene Symbol','name': 'Gene Name'}, inplace=True)
gene_info.to_csv('basic_gene_info.csv', index=False)

# Import gene tiers and create a dict
df2 = pd.read_csv('tiers_genes.csv')
data2 = df2
tiers_genes = pd.DataFrame(data2)
tier_data = tiers_genes['Tier']
frames = [gene_info, tier_data]  # The two existing dicts I want to merge
all_data = pd.concat(frames, axis=1)  # Merge the dicts together

# Add the all_data dict to a new .csv file using pandas
all_data.to_csv('gene_data.csv', index=False)
#
##------------------------------------------------------------------------------
## Use bokeh to create a visual for this data
##------------------------------------------------------------------------------
#"""
#"""
#
#all_data = all_data.copy()
#df = pd.read_csv('all_gene_data.csv', index_col=0)
#gene_data = df
#
## Normalize the data columns and sort
#gene_data.sort_values(by = 'Tier', inplace=True)
#
#column = []
#for x in gene_data.apply(tuple):
#    column.extend(x)
#
#
##n = all_data['Gene Name']
##t = all_data['Tier']
##s = all_data['Gene Symbol']
##summary = all_data['Gene Summary']
##all_data = all_data.copy()
##
##for name, tier, symbol in zip(n,t,s):
##    if tier == 1:
##        t1 = []
##        t1.append(symbol)
##        #print("Tier 1: "+  t1)
##    if tier == 3:
##        t3 = symbol
##        #print("Tier 3: "+  t3)
##    if tier == 2:
##        t2 = symbol
##        #print("Tier 2: "+  t2)
##    if tier == 0:
##        none_tiered = symbol
##        #print("None Tiered: "+  none_tiered)
#
##numbers = ["TIER 1", "TIER 2", "TIER 2", "NONE"]
#
#data = {
#  'tier': list(gene_data.Tier) * len(gene_data.columns),
#  'info':  [item for item in list(gene_data.columns) for i in range(len(gene_data.index))],
#  'column':   column,
#  'name': list(gene_data.index) * len(gene_data.columns)
#}
#
#output_file('test.html')
#hm = HeatMap(data, y='name', x='info',values='tier', title='Gene Viz', stat=None,
#             palette=RdYlGn4, legend=False, height=2500, spacing_ratio=0.7)
#
#
#
#
#hm.title.text_color = "black"
#hm.title.text_font = "times"
#hm.title.text_font_style = "bold"
##hm.add_layout(legend, 'right')
#hm.add_tools(HoverTool(tooltips=[("Gene Symbol", "@y"), ("Information Type", "@x"),("Tier", "@values")]))
#show(hm)


#n = all_data['Gene Name']
#t = all_data['Tier']
#s = all_data['Gene Symbol']
#summary = all_data['Gene Summary']
#all_data = all_data.copy()
#
#for name, tier, symbol in zip(n,t,s):
#    if tier == 1:
#        t1 = symbol
#        #print("Tier 1: "+  t1)
#    if tier == 3:
#        t3 = symbol
#        #print("Tier 3: "+  t3)
#    if tier == 2:
#        t2 = symbol
#        #print("Tier 2: "+  t2)
#    if tier == 0:
#        none_tiered = symbol
#        #print("None Tiered: "+  none_tiered)
#
#    numbers = ["TIER 1", "TIER 2", "TIER 2", "NONE"]
#
#    all_data = all_data.copy()
#    all_data["Entrez ID"] = all_data["Entrez ID"].astype(str)
#
#    tiers = [numbers[x-1] for x in zip(t1,t2,t3,none_tiered)]
#    group_range = [str(x) for x in range(1, 160)]
#
#
#    colormap = colormap = {
#        "Tier 1"         : "#a6cee3",
#        "Tier 2" : "#1f78b4",
#        "Tier 3"              : "#fdbf6f",
#        "None Tiered"                : "#b2df8a",}
#
#    source = ColumnDataSource(
#    data=dict(
#        tiers=[str(y) for y in zip(t1,t2,t3,none_tiered)],
#                symx=[str(x)+":0.1" for x in elements["group"]],
#        sym=all_data['Gene Symbol'],
#        name=all_data['Gene Name'],
#        summary=all_data['Gene Summary'],
#        entrez_id=all_data['Entrez ID'],
#        )
#    )
#
#    p = figure(title="Gene Information Map", tools="hover,save",
#               x_range=group_range, y_range=list(numbers))
#    p.plot_width = 1200
#    p.toolbar_location = None
#    p.outline_line_color = None
#
#    p.rect("tiers", 0.9, source=source,
#           fill_alpha=0.6, color="type_color")
#
#    text_props = {
#        "source": source,
#        "angle": 0,
#        "color": "black",
#        "text_align": "left",
#        "text_baseline": "middle"
#    }
#
#    p.text(x="symx", y="period", text="sym",
#           text_font_style="bold", text_font_size="15pt", **text_props)
#
#    p.text(x="symx", y="numbery", text="atomic_number",
#           text_font_size="9pt", **text_props)
#
#    p.text(x="symx", y="namey", text="name",
#           text_font_size="6pt", **text_props)
#
#    p.text(x="symx", y="massy", text="mass",
#           text_font_size="5pt", **text_props)
#
#    p.grid.grid_line_color = None
#
#    p.select_one(HoverTool).tooltips = [
#        ("name", "@name"),
#        ("atomic number", "@atomic_number"),
#        ("type", "@type"),
#        ("atomic mass", "@mass"),
#        ("CPK color", "$color[hex, swatch]:cpk"),
#        ("electronic configuration", "@electronic"),
#    ]
#
#    show(p)  # Change to save(p) to save but not show the HTML file
