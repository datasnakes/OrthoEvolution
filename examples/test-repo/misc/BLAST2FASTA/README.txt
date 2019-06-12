Description:

Blast2Fasta (blast2fasta.sh) is the encompassing of NCBI's BLAST to cds feature extraction.

1. Run GetGiLists.py to get a list of taxonomy ids for an Organisms.csv file. This also returns GI lists for each species.

2. Run AccCollect.py to run a the Standalone NCBI Blast client. It uses the GI lists from the prior step.

3. Run CleanCSV.py to clean the Master_Accession_File.csv with Blast.

4. Run xyz.py to create the gene genbank databases & extract cds fasta files.

5. Run GBK2Clustal.py to extract features & align them using Clustal Omega.

---------------------------------------------------------------------------------

Directions:

In order to run 'blast2fasta.sh', type 'sh blast2fasta.sh &' at the command line.

You may edit the bash script to include your email address so that the log file will be sent to you.

