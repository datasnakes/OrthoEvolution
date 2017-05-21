#!/usr/local/apps/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
File Name: blastntest.py
Description: This script inputs a list of organisms & genes. It tests a script.

@authors: Robert A. Gilmore & Shaurita D. Hutchins
Date Created: March 29, 2017
Project Name: Orthologs Project

Edited and updated for use with local/standalone NCBI BLAST 2.6.0
"""
import csv
import os
import subprocess
import shutil
import time  # Used to delay when dealing with NCBI server errors
from datetime import datetime as d
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline  # Used for Local Blasting.
from Manager.logit.logit import LogIt
from Orthologs.CompGenetics.ncbi_blast import BLASTAnalysis as BT


# TODO-ROB: Find packages for script timing and analysis
# TODO-ROB:  Add function for renaming and moving the builder files if the BLAST completed

class BLASTn(BT):
    def __init__(self, repo, user, project, research, research_type, template=None, save_data=True, **kwargs):
        """Inherit from the BLASTing Template."""
        super().__init__(template=template, repo=repo, user=user, project=project,
                         research=research, research_type=research_type, save_data=save_data, **kwargs)
        # # TODO-ROB Add taxon parameter
        # Manage Directories
        self.__home = Path(os.getcwd())
        self.__output_path = self.raw_data / Path('blast')  # Output directory
        self.__gi_list_path = self.__output_path / Path('gi_lists')
        self.__xml_path = self.__output_path / Path('xml')
        Path.mkdir(self.__output_path, parents=True, exist_ok=True)
        Path.mkdir(self.__gi_list_path, parents=True, exist_ok=True)
        Path.mkdir(self.__xml_path, parents=True, exist_ok=True)
        # # Initialize Logging
        # self.__blastn_log = LogIt.blastn()
        #df = LogIt()
        #self.date_format = df.date_format
        # self.get_time = time.time  # To get the time use 'get_time()'
        # TODO-ROB:  Add a query organism variable
        self.query_gi_dict = {}
        self.removed_genes = []
        # TODO-ROB:  Set up blast config logger, blasting logger, and post blast analysis logger
        # ------------------------------------------------------------------------------
        self.blastn_log.info("These are the organisms: " + str(self.org_list))
        self.blastn_log.info("These are the genes: " + str(self.gene_list))
        self.blastn_log.info("These are the taxonomy ids: \n\n\n" + str(self.taxon_ids))
        # ---------------------------------------------------------------------
        # For completed blast files
        self.complete_file = self.project + '_MAF.csv'
        self.complete_file_path = self.data / Path(self.complete_file)
        self.complete_time_file = self.project + '_TIME.csv'
        self.complete_time_file_path = self.data / Path(self.complete_time_file)

    @staticmethod
    def map_func(hit):
        """The map function for formatting hit id's.
        This will be used later in the script."""

        hit.id1 = hit.id.split('|')[3]
        hit.id2 = hit.id.split('|')[1]
        hit.id = hit.id[:-2]
        return hit

    def blast_config(self, query_align, query_organism, auto_start=False):
        """This function configures everything for a BLAST.
        First the accession file, and gene list is configured."""
        # os.chdir(str(self.__output_path))
        self.blastn_log.info('***********************************BLAST CONFIG START***********************************')
        self.blastn_log.info('***********************************BLAST CONFIG START***********************************')
        self.blastn_log.info('***********************************BLAST CONFIG START***********************************\n\n\n')

        self.blastn_log.info('Configuring the accession file...')
        # Update the gene_list based on the existence of a incomplete blast file
        gene_list = self.blast_file_config(self.building_file_path)
        if gene_list is not None:
            start = len(self.blast_human) - len(gene_list)  # Number of genes already BLASTed
            query_align = self.blast_human[start:]  # Reconfigure query_align to reflect the existing accession info
            new_gene_list = gene_list  # Reconfigure the gene_list to reflect the existing accession info
        else:
            new_gene_list = self.gene_list

        # Create GI lists
        self.blastn_log.info("Configuring GI list using the taxonomy id and the blastdbcmd tool.")
        #self.gi_list_config()
        # Get GI (stdout) and query sequence (FASTA format)
        self.blastn_log.info("Generating directories.")
        self.blastn_log.info("Extracting query gi number to stdout and "
                             "query refseq sequence to a temp.fasta file from BLAST database.")
        # Iterate the query accessions numbers
        for query in query_align:
            # os.chdir(str(self.__output_path))
            gene = self.acc_dict[query][0]
            gene_path = self.__xml_path / Path(gene)
            org = self.acc_dict[query][1]
            # Create the proper directories for each gene
            try:
                Path.mkdir(gene_path)
                self.blastn_log.info("Directory Created: %s" % gene)
                self.blastn_log.info("\n")
                #os.chdir(gene)
            except FileExistsError:
                self.blastn_log.info("Directory already exists: %s" % gene)
                #os.chdir(gene)


            # Save sequence data in FASTA file format and print the gi number to stdout with a custom BLAST extraction
            # https://www.ncbi.nlm.nih.gov/books/NBK279689/#_cookbook_Custom_data_extraction_and_form_
            # TODO-ROB:  TODO-SHAE:Combine these BLAST extractions???
            fmt = {'query': query, 'temp fasta': str(gene_path / Path('temp.fasta'))}
            fasta_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %f -out {temp fasta}".format(**fmt)
            fasta_status = subprocess.call([fasta_setup], shell=True)
            gi_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %g".format(**fmt)
            gi_status = subprocess.call([gi_setup], shell=True)
            # TODO-ROB:  Add function to add the gi numbers to the dataframe/csv-file, and add a check function to see if thats already there
            # Check the status of the custom blast data extraction
            if gi_status == 0 or fasta_status == 0:  # Command was successful.
                if gi_status != 0:
                    self.blastn_log.error("GI number for %s not found in the BLAST extraction" % query)  # Log it.
                    # TODO-ROB: TODO-SHAE: Is this the correct move below???
                    self.blastn_log.error("Removing %s from the BLAST list..." % gene)
                    self.gene_list.remove(gene)
                    self.removed_genes.append(gene)
                    continue
                if fasta_status != 0:
                    self.blastn_log.error("FASTA sequence for %s not found in the BLAST extraction" % query)
                    self.blastn_log.error("Removing %s from the BLAST list..." % gene)
                    self.gene_list.remove(gene)
                    self.removed_genes.append(gene)
                    continue
                pass
            else:
                self.blastn_log.error("FASTA sequence and GI number for %s not found in the custom BLAST extraction." % query)
                self.blastn_log.error("Removing %s from the BLAST list..." % gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)
                continue

            # Get the gi number from stdout, format it, and add it to the gi dictionary
            gi = subprocess.check_output([gi_setup], shell=True)
            gi = gi.strip()
            gi = gi.decode('utf-8')
            gi = str(gi)
            gi = gi.replace("'", "")
            self.query_gi_dict[gene] = gi

        new_query_align = query_align
        self.blastn_log.info('Configured query accession list: %s' % new_query_align)
        self.blastn_log.info('Configured gene list: %s\n\n\n' % new_gene_list)
        self.blastn_log.info('************************************BLAST CONFIG END************************************')
        self.blastn_log.info('************************************BLAST CONFIG END************************************')
        self.blastn_log.info('************************************BLAST CONFIG END************************************\n\n\n')
        if auto_start is True:
            # Automatically begin BLASTING after the configuration
            self.blasting(genes=new_gene_list, query_organism=query_organism, pre_configured=auto_start)
        else:
            # Manually begin BLASTING and return the new gene and new query lists
            return new_gene_list

    def gi_list_config(self):
        # TODO-ROB THis is for development / testing
        """This script is designed to create a gi list based on the refseq_rna database
        for each taxonomy id on the MCSR. It will also convert the gi list into a
        binary file which is more efficient to use with NCBI's Standalone Blast tools."""
        print('gi_list_config')
        cd = os.getcwd()
        os.chdir(str(self.__gi_list_path))
        taxids = self.taxon_ids
        pd.Series(taxids).to_csv('taxids.csv', index=False)
        #subprocess.call(['qsub %s' % str(self.__gi_list_path / Path('get_gi_lists.pbs'))], shell=True)
        p = subprocess.Popen(['qsub %s' % str(self.__gi_list_path / Path('get_gi_lists.pbs'))], shell=True)
        p.wait()
        print('Done submitting jobs')
        gi_flag = True
        while gi_flag == True:
            try:
                subprocess.check_output(['pidof', 'getgilists'])
                gi_flag = False
            except subprocess.CalledProcessError:
                gi_flag = True
                print('Waiting for the gi BLAST to finish....')
                time.sleep(30)
        print('Multiprocessing complete')
        os.chdir(cd)
    #     with Pool(processes=20) as p:
    #         cd = os.getcwd()
    #         os.chdir(self.__gi_list_path)
    #         p.map(self.gi_split, taxids)
    #         os.chdir(cd)
    #         self.blastn_log.inf("The GI lists have been created.")
    #
    # @staticmethod
    # def gi_split(ID):
    #     """ This function uses the blastdbcmd tool to get gi lists. It then uses the
    #     blastdb_aliastool to turn the list into a binary file."""
    #     # Create dictionary for formatting
    #     ID = str(ID)
    #     gi_text_file = "%sgi.txt" % ID
    #     gi_binary_file = "%sgi" % ID
    #     fmt = {'id': ID, 'text file': gi_text_file, 'binary file': gi_binary_file}
    #     # TODO untested
    #     print("Current taxonomy ID: %s" % ID)
    #     # Use the accession #'s and the blastdbcmd tool to generate gi lists based on Organisms/Taxonomy ID's.
    #     os.system("blastdbcmd -db refseq_rna -entry all -outfmt '%g %T' | awk ' {{ if ($2 == {id}) {{ print $1 }} }} ' > {text file}".format(**fmt))
    #     print("Text File generated")
    #     # Convert the .txt file to a binary file using the blastdb_aliastool.
    #     os.system("blastdb_aliastool -gi_file_in {text file} -gi_file_out {binary file}".format(**fmt))
    #     # Remove the gi.text file
    #     print("Binary file generated")
    #     os.system("rm -r {text file}".format(**fmt))
    #     print("Removing text file")

    def blast_file_config(self, file):
        """This function configures different files for new BLASTS.
        It also helps recognize whether or not a BLAST was terminated
        in the middle of the dataset.  This removes the last line of
        the accession file if it is incomplete."""

        output_dir_list = os.listdir(self.__output_path)  # Make a list of files
        # If the file exists then make a gene list that picks up from the last BLAST
        if file in output_dir_list:
            with open(file, 'r') as fi:
                f = csv.reader(fi)
                count = - 1
                for row in f:
                    count += 1
                    ending = row
                gene = ending[1]
                taxid = self.taxon_ids[count]
                org = self.org_list[(len(ending) - 2)]

                ncbi = str("""result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna", strand="plus", 
                evalue=0.001, out="%s_%s.xml", outfmt=5, gilist=%s + "gi", max_target_seqs=10, task="blastn")"""
                           % (gene, org, taxid))
                self.blastn_log.warning("An incomplete accession file was produced from the previous BLAST,"
                                        "which was terminated midway through the procedure.")

                self.blastn_log.info("The last row looks like: \n\t%s\n\t%s\n" % (self.header, ending))
                self.blastn_log.info("The BLAST ended on the following query: \n%s" % ncbi)
                if len(ending) < len(self.header):
                    self.blastn_log.info("Restarting the BLAST for the previous gene...")
                    count = count - 1
                # The continued gene list starts with the previous gene.
                continued_gene_list = list(x for i, x in enumerate(self.gene_list, 1) if i > count)
            return continued_gene_list
        # If the file doesn't exist return nothing
        else:
            self.blastn_log.info("A new BLAST started at %s" % self.get_time())
            return None

    def blast_xml_parse(self, xml_file, gene, organism):
        self.blastn_log.info("Parsing %s to find the best accession number." % xml_file)
        # Parse the XML file created by the BLAST
        maximum = 0
        with open(xml_file, 'r') as blast_xml:
            blast_qresult = SearchIO.read(blast_xml, 'blast-xml')
            mapped_qresult = blast_qresult.hit_map(self.map_func)
            for hit in mapped_qresult:
                for hsp in hit.hsps:
                    # Find the highest scoring hit for each gene
                    if hsp.bitscore_raw > maximum:
                        # If the gene is a predicted non-coding RefSeq gene then go the the next hit
                        # https://en.wikipedia.org/wiki/RefSeq
                        # TODO-ROB:  TODO-SHAE:  Add another check here???
                        if "xr" in str(hit.id.lower()):
                            self.blastn_log.info("Encountered a predicted(X*_00000) "
                                                 "non-coding (XR_000000)(ncRNA) RefSeq gene.  Moving to the next hit.")
                            break
                        else:
                            # If the gene is acceptable then add it to the gene list
                            # Lowercase means LOC1223214 is the name of the gene
                            # TODO-ROB:  Change this check
                            maximum = hsp.bitscore_raw
                            if gene.lower() in hit.description.lower():
                                accession = hit.id1
                                gi = hit.id2
                                raw_bitscore = hsp.bitscore_raw
                                description = hit.description
                            else:
                                accession = hit.id1.lower()
                                gi = hit.id2
                                raw_bitscore = hsp.bitscore_raw
                                description = hit.description
            self.blastn_log.info("Finished parsing the file.")
            self.blastn_log.info("The best accession has been selected from the BLAST xml record.")
            self.blastn_log.info("Accession:  %s" % accession)
            self.blastn_log.info("GI number: %s" % gi)
            self.blastn_log.info("Raw bitscore: %s" % raw_bitscore)
            self.blastn_log.info("Description: %s" % description)
            self.add_accession(gene, organism, accession)

# *********************************************************************************************************************
# *********************************************************************************************************************
# *********************************************************************************************************************

    def blasting(self, genes=None, query_organism=None, pre_configured=False):
        # Configure the BLAST
        if pre_configured is False:
            query = self.df[query_organism].tolist()
            genes = self.blast_config(query_align=query, query_organism=query_organism, auto_start=True)
        elif pre_configured is True:
            genes = genes

        self.blastn_log.info("------------------------------------------------------------------")
        self.blastn_log.info("The script name is %s" % os.path.basename(__file__))
        self.blastn_log.info("The script began on %s" % str(d.now().strftime(self.date_format)))
        self.blastn_log.info("------------------------------------------------------------------")

        # Steps to a bulk blast
        # 1.  Iterate the gene of interest
        # 2.  Iterate the organisms of interest
        # 3.  BLAST
        for gene in genes:
            for organism in self.org_list:
                # Skip the query organism
                if organism == query_organism:
                    continue
                # Initialize output variables
                gene_path = self.__xml_path / Path(gene)
                files = os.listdir(gene_path)
                xml = '%s_%s.xml' % (gene, organism)
                xml_path = gene_path / Path(xml)

                # Initialize configuration variables
                taxon_id = self.taxon_dict[organism]
                taxon_gi_file = str(taxon_id + "gi")
                taxon_gi_path = self.__gi_list_path / Path(taxon_gi_file)
                taxgi_dest_path = gene_path / Path(taxon_gi_file)
                if xml in files:
                    self.blast_xml_parse(xml_path, gene, organism)
                else:
                    self.blastn_log.warning("\n\n\n*******************BLAST START*******************")
                    start_time = self.get_time()
                    self.blastn_log.info("The start time is %s" % start_time)
                    self.blastn_log.info("The current gene is %s (%s)." % (gene, self.tier_dict[gene]))
                    self.blastn_log.info("The current organisms is %s (%s)." % (organism, taxon_id))

                    with open(xml_path, 'w') as blast_xml:
                        # Create a copy of the gi list file per taxonomy id to be used in blast
                        fmt = {'src': str(taxon_gi_path), 'dst': str(taxgi_dest_path)}
                        os.system("cp {src} {dst}".format(**fmt))
                        # Set up blast parameters
                        taxgi_dest_path = str(taxgi_dest_path)
                        query_seq_path = str(gene_path / Path('temp.fasta'))
                        # Use Biopython's NCBIBlastnCommandline tool
                        result_handle1 = NcbiblastnCommandline(query=query_seq_path, db="refseq_rna",
                                                               strand="plus", evalue=0.001,  # DONT GO LOWER
                                                               outfmt=5, gilist=taxgi_dest_path,
                                                               max_target_seqs=10, task="blastn")
                        # Capture blast data
                        stdout_str, stderr_str = result_handle1()
                        blast_xml.write(stdout_str)
                        end_time = self.get_time()
                        elapsed_time = end_time - start_time
                        self.blastn_log.info("%s was create." % blast_xml.name)
                        self.blastn_log.info("The end time is %s." % end_time)
                        self.blastn_log.info("The BLAST took %s." % elapsed_time)
                    self.blastn_log.warning("********************BLAST END********************\n\n\n")
                    self.add_blast_time(gene, organism, start_time, end_time)
                    self.blast_xml_parse(xml, gene, organism)
                # If the BLAST has gone through all orthologs then create a master accession file.
                if gene == genes[-1] and organism == self.org_list[-1]:
                    shutil.copy(str(self.building_file_path), str(self.complete_file_path))
                    shutil.copy(str(self.building_time_file_path), str(self.complete_time_file_path))
                    if self.save_data is True:
                        self.post_blast_analysis(self.project)
                    # TODO-ROB Archive function here
