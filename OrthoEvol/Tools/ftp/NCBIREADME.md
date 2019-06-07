#### Contents of the /blast/db/ directory

The names of these databases and their contents are listed below.

 File Name                    | Content Description
|-----------------------------|------------------------------------------------|
16SMicrobial.tar.gz           | Bacterial and Archaeal 16S rRNA sequences from BioProjects 33175 and 33117
FASTA/                        | Subdirectory for FASTA formatted sequences
README                        | README for this subdirectory (this file)
v5/                           | Subdirectory for preformatted version 5 blast databases
Representative_Genomes.*tar.gz| Representative bacterial/archaeal genomes database
cdd_delta.tar.gz              | Conserved Domain Database sequences for use with stand alone deltablast
cloud/	                      | Subdirectory of databases for BLAST AMI; see http://1.usa.gov/TJAnEt
env_nr.*tar.gz                | Protein sequences for metagenomes
env_nt.*tar.gz                | Nucleotide sequences for metagenomes
est.tar.gz                    | This file requires est_human.*.tar.gz, est_mouse.*.tar.gz, and est_others.*.tar.gz files to function. It contains the est.nal alias so that searches against est (-db est) will include est_human, est_mouse and est_others.
est_human.*.tar.gz            | Human subset of the est database from the est division of GenBank, EMBL and DDBJ.
est_mouse.*.tar.gz            | Mouse subset of the est database
est_others.*.tar.gz           | Non-human and non-mouse subset of the est database
gss.*tar.gz                   | Sequences from the GSS division of GenBank, EMBL, and DDBJ
htgs.*tar.gz                  | Sequences from the HTG division of GenBank, EMBL, and DDBJ
human_genomic.*tar.gz         | Human RefSeq (NC_######) chromosome records with gap adjusted concatenated NT_ contigs
nr.*tar.gz                    | Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq
nt.*tar.gz                    | Partially non-redundant nucleotide sequences from all traditional divisions of GenBank, EMBL, and DDBJ excluding GSS, STS, PAT, EST, HTG, and WGS.
other_genomic.*tar.gz         | RefSeq chromosome records (NC_######) for non-human organisms
pataa.*tar.gz                 | Patent protein sequences
patnt.*tar.gz                 | Patent nucleotide sequences. Both patent databases are directly from the USPTO, or from the EPO/JPO via EMBL/DDBJ
pdbaa.*tar.gz                 | Sequences for the protein structure from the Protein Data Bank
pdbnt.*tar.gz                 | Sequences for the nucleotide structure from the Protein Data Bank. They are NOT the protein coding sequences for the corresponding pdbaa entries.
refseq_genomic.*tar.gz        | NCBI genomic reference sequences
refseq_protein.*tar.gz        | NCBI protein reference sequences
refseq_rna.*tar.gz            | NCBI Transcript reference sequences
sts.*tar.gz                   | Sequences from the STS division of GenBank, EMBL, and DDBJ
swissprot.tar.gz              | Swiss-Prot sequence database (last major update)
taxdb.tar.gz                  | Additional taxonomy information for the databases listed here providing common and scientific names
tsa_nt.*tar.gz                | Sequences from the TSA division of GenBank, EMBL, and DDBJ
vector.tar.gz                 | Vector sequences from 2010, see Note 2 in section 4.

#### Contents of the /blast/db/v5/ directory

 File Name                    | Content Description
|-----------------------------|------------------------------------------------|
blastdbv5.pdf                 | A pdf of examples for blastdbv5.
nr_v5.*tar.gz                 | Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq
nt_v5.*tar.gz                 | Partially non-redundant nucleotide sequences from all traditional divisions of GenBank, EMBL, and DDBJ excluding GSS, STS, PAT, EST, HTG, and WGS.
refseq_rna_v5.*tar.gz         | NCBI Transcript reference sequences
swissprot_v5.tar.gz           | Swiss-Prot sequence database (last major update)
taxdb.tar.gz                  | Additional taxonomy information for the databases listed here providing common and scientific names
