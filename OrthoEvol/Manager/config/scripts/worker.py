import os
import sys
from pathlib import Path

from Archive.filter import FilteredAlignment
from OrthoEvol.Orthologs.Phylogenetics.IQTree.best_tree import FilteredTree
from OrthoEvol.Orthologs.Phylogenetics.PAML.codeml import CodemlRun

raw_data_path = sys.argv[2]  # path/raw_data/gene/
os.chdir(raw_data_path)

index_path = str(Path(os.getcwd()).parent / Path('index'))  # path/index
gene = sys.argv[1]  # gene

paml_control_file = str('paml.ctl')  # path/index/paml_control_file

iqtree_newick = str(gene + '_iqtree.nwk')  # path/raw_data/gene/iqtree_newick

# 1.  Run a Guidance2 based alignment filter on the nucleotide sequences for at least 2 iterations.
na_fasta = str(gene + '.ffn')  # path/raw_data/gene/na_fasta
# 2.  Filter the amino acid sequences based on the first step.
aa_fasta = str(gene + 'faa')  # path/raw_data/gene/aa_fasta
# 3.  Run a Guidance2 based alignment filter on the filtered amino acid sequences.
# 4.  Using Pal2Nal create a FASTA formatted alignment for IQTree.
p2n_iqtree_aln = str(gene + '_P2N_na.iqtree.aln')  # path/raw_data/gene/p2n_iqtree_aln
# 5.  Using Pal2Nal create a PAML formatted alignment for PAML.
p2n_paml_aln = str(gene + '_P2N_na.paml.aln')  # path/raw_data/gene/p2n_paml_aln
qca = FilteredAlignment(na_fasta, aa_fasta, na_colCutoff=0.1)

# 6. Run IQ-Tree on the codon alignment to get the best tree.
qct = FilteredTree(p2n_iqtree_aln)

# 7. Run Codeml to complete the phylogenetic analysis
x = CodemlRun(p2n_paml_aln, iqtree_newick, paml_control_file)
# 8. Run ggtree to generate visualizations
x.cml.run(verbose=True)
