Window Masker Setup
================

Creating a masked BLAST database
--------------------------------

Creating a masked BLAST database is a two step process:

1.  Generate the masking data using a sequence filtering utility like windowmasker or dustmasker

2.  Generate the actual BLAST database using makeblastdb

### Create a counts file

``` bash
windowmasker -in refseq_rna -infmt blastdb -mk_counts -parse_seqids -sformat obinary -out refseq_rna.counts
```

With the counts file we can then proceed to create the file containing the masking information as follows:

``` bash
windowmasker -in refseq_rna -infmt blastdb -ustat refseq_rna.counts -outfmt maskinfo_asn1_bin -parse_seqids -out refseq_rna.asnb
```

### Create BLAST database with the masking information

Using the masking information data files generated in the previous steps, we can create BLAST database with masking information incorporated.

Note: we should use â€œ-parse\_seqidsâ€ in a consistent manner â€“ either use it in both steps or not use it at all.

``` bash
makeblastdb -in refseq_rna -input_type blastdb -dbtype nucl -parse_seqids -mask_data refseq_rna.asnb -out refseq_rna -title "Refseq RNA Masked Database"
```

``` bash
blastdbcmd -db refseq_rna -info
```
