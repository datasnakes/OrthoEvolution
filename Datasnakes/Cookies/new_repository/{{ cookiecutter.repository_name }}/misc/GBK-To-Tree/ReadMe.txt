1. Extracts and saves a feature from a genbank file.

2. Translate the feature (sequence) to the amino acid sequence.

3. Uses that file (fasta formatted) with Clustal Omega to produce an alignment (in phylip format).

4. Executes the PHYLIP program to produce maximum likelihood trees and distance matrices.

5. Uses the output of the prior program (phylip format) with PAML to generate analysis of sequences.