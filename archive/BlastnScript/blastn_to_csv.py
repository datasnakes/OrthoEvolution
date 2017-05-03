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
# Modules as custom
import logging as log
# Modules used
import os
import subprocess
import time  # Used to delay when dealing with NCBI server errors
from datetime import datetime as d

import pandas as pd  # Used for dealing with CSV files
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline  # Used for Local Blasting.
# Proprietary Modules
from manager.ortho_analysis import OrthologAnalysis

from Orthologs.manager.blast_analysis import BLASTAnalysis as BT

# TODO-ROB: Find packages for script timing and analysis

# Set up the blastn logger & log file
format1 = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
format2 = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
LOGFORMAT = '%(name)s - [%(levelname)-2s]: %(message)s'
log.basicConfig(level=log.DEBUG,
                format=LOGFORMAT,
                filename="logs/accessions2blastxml_%s.log" % str(d.now().strftime(format2)))
blastn_log = log.getLogger('Blastn')

# Write basic information to the log
blastn_log.info("------------------------------------------------------------------")
blastn_log.info("The script name is %s" % os.path.basename(__file__))
blastn_log.info("The script began on %s" % str(d.now().strftime(format1)))
blastn_log.info("------------------------------------------------------------------")

# Initializations
check = 1  # Variable used to check the output
start_time = time.time()  # Variable used to check the processing time
get_time = time.time  # To get the time use 'get_time()'
# Initialize from the template master accession file and taxon file
template = 'MAFV3.1.csv'
data = BT(template=template)
"""The template accession file contains the following information...
Headers: Tier, Gene, Org_1, ... , Org_n
Rows:  tier_name, gene_name, Accession_1, ... , Accession_n"""
header = data.header
org_list = data.org_list
gene_list = data.gene_list
taxon_ids = data.taxon_ids
"""The taxon ids file contains one id on each line"""

# Manage directories
h = os.getcwd()  # Home directory
output = 'data/blastn/blast-xml-output/'  # Output directory
processed = 'data/processed/'  # Processed data directory
os.makedirs(output, exist_ok=True)  # only make blast-xml-output dir since blastn dir exists

# ------------------------------------------------------------------------------
blastn_log.info("These are the organisms: " + str(org_list))
blastn_log.info("These are the genes: " + str(gene_list))
blastn_log.info("These are the taxonomy ids: " + str(taxon_ids))
# ------------------------------------------------------------------------------
os.chdir(output)  # Change or move to the output directory
output_dir_list = os.listdir()  # Make a list of files

# ------------------------------------------------------------------------------
# Check to see if the master ACCESSION file is in the home directory,
# and then either add the header or count the number of rows that already
# exist in order to skip to the most recently called gene.


def blast_file_config(file):
    """This function configures different files for new BLASTS.
    It also helps recognize whether or not a BLAST was terminated
    in the middle of the dataset.  This removes the last line of
    the accession file if it is incomplete."""
    if file in output_dir_list:
        with open(file, 'r') as fi:
            f = csv.reader(fi)
            count = -1
            for row in f:
                count += 1
                ending = row
            gene = ending[1]
            taxid = taxon_ids[count]
            org = org_list[(len(ending) - 2)]

            ncbi = str("""result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna", strand="plus",
            evalue=0.001, out="%s_%s.xml", outfmt=5, gilist=%s + "gi", max_target_seqs=10, task="blastn")"""
                       % (gene, org, taxid))
            blastn_log.warning("An incomplete accession file was produced from the previous BLAST,"
                               "which was terminated midway through the procedure.")

            blastn_log.info("The last row looks like: \n\t%s\n\t%s\n" % (header, ending))
            blastn_log.info("The BLAST ended on the following query: \n%s" % ncbi)
            if len(ending) < len(header):
                blastn_log.info("Removing the last line...")
                count = count -1
            continued_org_list = list(x for i, x in enumerate(gene_list, 1) if i > count)
            return continued_org_list
    else:
        blastn_log.info("A new BLAST started at %s" % get_time())

if 'Master_Accession_File.csv' in output_dir_list:
    MAF = open('Master_Accession_File.csv')
    MAF = csv.reader(MAF)
    row_count_A = sum(1 for row in MAF) - 1
    blastn_log.info("row_count_A: " + str(row_count_A))

else:
    MAF = open('Master_Accession_File.csv', 'w', newline='')
    org_row = csv.writer(MAF, dialect='excel')
    org_row.writerow(['Tier'] + ['Gene'] + org_list)
    MAF.close()
    row_count_A = 0
    blastn_log.info("2_A")

#
# DEPRECATED------------------------------------------------------------------------------
# Check to see if the master GI file is in the home directory,
# and then either add the header or count the number of rows that already
# # exist in order to skip to the most recently called gene
#
# if 'Master_GI_File.csv' in output_dir_list:
#     MGF = open('Master_GI_File.csv')
#     MGF = csv.reader(MGF)
#     row_count_G = sum(1 for row in MGF) - 1
#     blastn_log.info("row_count_G: " + str(row_count_G))
#     # input("Is this an ok row_count for the GI File?")
# else:
#     MGF = open('Master_GI_File.csv', 'w', newline = '')
#     org_row = csv.writer(MGF, dialect='excel')
#     org_row.writerow(['Tier'] + ['Gene'] + org_list)
#     MGF.close()
#     row_count_G = 0
#     blastn_log.info("2_G")

# DEPRECATED------------------------------------------------------------------------------
# Check to see if the Time_Record_ExistingFiles file is in the home directory,
# and then either add the header or count the number of rows that already exist
# in order to skip to the most recently called gene

if 'Time_Record_ExistingFiles.csv' in output_dir_list:
    TREF = open('Time_Record_ExistingFiles.csv')
    TREF = csv.reader(TREF)
    row_count_TREF = sum(1 for row in TREF) - 1
    blastn_log.info("row_count_TREF: " + str(row_count_TREF))
    # input("Is this an ok row_count for the Time_Record_ExistingFiles File?")
else:
    TREF = open('Time_Record_ExistingFiles.csv', 'w', newline='')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(['Tier'] + ['Gene'] + org_list)
    TREF.close()
    row_count_TREF = 0
    blastn_log.info("2_TREF")

# DEPRECATED------------------------------------------------------------------------------
# Check to see if the Time_Record_BLASTingFiles file is in the home directory,
# and then either add the header or count the number of rows that already exist
# in order to skip to the most recently called gene

if 'Time_Record_BLASTingFiles.csv' in output_dir_list:
    TRBF = open('Time_Record_BLASTingFiles.csv')
    TRBF = csv.reader(TRBF)
    row_count_TRBF = sum(1 for row in TRBF) - 1
    blastn_log.info("row_count_TRBF: " + str(row_count_TRBF))
else:
    TREF = open('Time_Record_BLASTingFiles.csv', 'w', newline='')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(['Tier'] + ['Gene'] + org_list)
    TREF.close()
    row_count_TREF = 0
    blastn_log.info("2_TRBF")

# DEPRECATED------------------------------------------------------------------------------
# For writing Accessions/GIs to the file
gene_list_A = []
gene_list_B = []
file_count = 0
BA_count = 0
GI_count = 0

# ------------------------------------------------------------------------------
# Define the map function for hit id's.
# This will be used later in the script.
def map_func(hit):
    hit.id1 = hit.id.split('|')[3]
    hit.id2 = hit.id.split('|')[1]
    hit.id = hit.id[:-2]
    return hit

# ------------------------------------------------------------------------------
# BLAST the Homo Sapiens Accession number & filter results via the Organism
# The 1st 'for' loop parses through the individual genes in our data files
os.chdir(h)  # Change to home directory

# 1st column - tiers, 2nd column - genes, 3rd column - Human accession numbers
# TODO-ROB: CompGenAnalysis here
file1 = csv.reader(open('data/initial-data/homo_sapiens_accessions.csv'))
os.chdir(output)

Acc_count = 0  # The Accession count. Lists start at 0.

for Accession in file1:
    Acc_count = Acc_count + 1

#  Begin listing I/O information
# ------------------------------------------------------------------------------
    # TODO-ROB:  CompGenAnalysis here
    blastn_log.info('#' + (50 * '-'))
    blastn_log.info('The following contains information about the BLAST input and output:')
    blastn_log.info('Human Accession: %s' % Accession[2])
    blastn_log.info("Gene of Interest: %s" % Accession[1])
    blastn_log.info("Gene Tier: %s" % Accession[0])

    blastn_log.info("Current Directory: " + str(os.getcwd()))  # Check the Python Shell for proper directory

    unfound_genes = []
    unfound_accs = []
    unfound_tier = []

# DEPRECATED                    DEPRECATED
    #     gene_list_A = []
    #     gene_list_G = []
    #     gene_list_TREF = []
    #     gene_list_TRBF = []

    # # Append gene tiers to the row
    # gene_list_A.append(str(Accession[0]))
    # gene_list_G.append(str(Accession[0]))
    # gene_list_TREF.append(str(Accession[0]))
    # gene_list_TRBF.append(str(Accession[0]))
    #
    # # Append gene names to the row
    # gene_list_A.append(str(Accession[1]))
    # gene_list_G.append(str(Accession[1]))
    # gene_list_TREF.append(str(Accession[1]))
    # gene_list_TRBF.append(str(Accession[1]))
# DEPRECATED                        DEPRECATED
# Create a folder (for target genes) unless it already exists.
# If it does exist, then change to that folder to see
# what XML files (for target gene/organisms) are present.

    x = str(Accession[1]) + "_BlastXML"  # Used for making/changing directories

    try:
        os.mkdir('%s' % (x))
        blastn_log.info("Directory Created: %s" % (x))
        os.chdir(x)
        blastn_log.info("Current Directory: " + str(os.getcwd()))
        blastn_log.info("\n")

    except FileExistsError:
        blastn_log.info("Directory already exists: %s" % (x))
        os.chdir(x)
        blastn_log.info("Current Directory: " + str(os.getcwd()))

    blastn_log.info(os.getcwd())  # Check the Python Shell for proper output
    blastn_log.info('#' + (50 * '-') + '\n' + '\n')

#  End listing I/O information

# ------------------------------------------------------------------------------
# The 2nd 'for' loop is for parsing through the individual organisms(which correspond
# to taxonomy id's) in our data files, so that very specific BLAST data can be
# retrieved and stored in the above directory.
    for Organism, TID in zip(org_list, taxon_ids):
        file_count = file_count + 1
        maximum = 0

        if Organism == 'Homo_sapiens':
            gene_list_A.append(str(Accession[2]))
            BA_count = BA_count + 1
            cmd = "blastdbcmd -entry " + str(Accession[2]) + " -db refseq_rna -outfmt %f -out temp.fasta"
            getgi = "blastdbcmd -entry " + str(Accession[2]) + " -db refseq_rna -outfmt %g"
            cmdstatus = subprocess.call([getgi], shell=True)

            if cmdstatus == 0:  # Command was successful.
                pass  # Continue through script.
            else:
                unfound_genes.append(str(Accession[1]))  # add unfound gene to list
                unfound_accs.append(str(Accession[2]))
                unfound_tier.append(str(Accession[0]))
                blastn_log.error("Entry not found: %s" % str(Accession[2]))  # Log it.

                # Remove the entry/gene/accession not found in the db from the lists.
                gene_list_A.remove(str(Accession[1]))
                gene_list_G.remove(str(Accession[1]))
                gene_list_TREF.remove(str(Accession[1]))
                gene_list_TRBF.remove(str(Accession[1]))
                break

            gi = subprocess.check_output([getgi], shell=True)

            gi = gi.strip()
            gi = gi.decode('utf-8')
            gi = str(gi)
            gi = gi.replace("'", "")
            gene_list_G.append(gi)
            GI_count = GI_count + 1

            blastn_log.info("----%s seconds----" % (time.time() - start_time))
            timer = str(time.time() - start_time)
            gene_list_TREF.append('')
            gene_list_TRBF.append(timer)
            pass

# ------------------------------------------------------------------------------
# Skip the first 3 items in the org_list and update the Best Accession Count and Best GI Count
        else:  # or Organism == 'Homo sapiens' or  Organism == 'Macaca mulatta':
            if Organism == '':
                BA_count = file_count
                GI_count = file_count
                continue
            xmlfile_list = os.listdir()
            s = "%s_%s.xml" % (Accession[1], Organism) #TODO-ROB REMOVE????

# ------------------------------------------------------------------------------
# If a gene/organism XML file has already been established, parse through it.
            if s in xmlfile_list:
                blastn_log.info(Organism)
                blastn_log.info(org_list)
                with open("%s_%s.xml" % (Accession[1], Organism)) as result_handle2:
                    blast_qresult = SearchIO.read(result_handle2, 'blast-xml')
                    mapped_qresult = blast_qresult.hit_map(map_func)
                    for hit in mapped_qresult:
                        for hsp in hit.hsps:
                            blastn_log.info(hit.id1)
                            blastn_log.info(hit.id2)
                            blastn_log.info(hsp.bitscore_raw)
                            blastn_log.info(hit.description + '\n')
                            if hsp.bitscore_raw > maximum:
                                if "xr" in str(hit.id.lower()):
                                    blastn_log.info("Encountered a pseudogene.  Moving to the next hit.")
                                    break
                                else:
                                    maximum = hsp.bitscore_raw
                                    if x.lower() in hit.description.lower():
                                        Best_Accession = hit.id1
                                        Best_GI = hit.id2
                                    else:
                                        Best_Accession = hit.id1.lower()
                                        Best_GI = hit.id2
                                    blastn_log.info(Best_Accession + ' ' + Best_GI + " has the highest bitscore!!!!")
                                    blastn_log.info(hsp.hit)
                                    gene_list_A.append(Best_Accession)
                                    gene_list_G.append(Best_GI)
                                    BA_count = BA_count + 1
                                    GI_count = GI_count + 1
                    blastn_log.info("BA_count = %s" % BA_count)
                    blastn_log.info("GI_count = %s" % GI_count)
                    blastn_log.info("file_count = %s" % file_count)
                    if file_count != BA_count:
                        gene_list_A.append('')
                        BA_count = file_count
                    if file_count != GI_count:
                        gene_list_G.append('')
                        GI_count = file_count
                    blastn_log.info("----%s seconds----" % (time.time() - start_time))
                    timer = str(time.time() - start_time)
                    gene_list_TREF.append(timer)
                    gene_list_TRBF.append('')

                blastn_log.info(s + " already exists.  Next organism.")
                continue

# ------------------------------------------------------------------------------
# If the gene/organism xml file hasn't already been established,
# then BLAST the Homo sapiens accession number,
# and store the blast report in an XML file.
            else:
                # TODO-ROB:  There is an error causing the csv file to write an extra space
                # TODO-ROB:  The 'if' is fine.  LINE206 in the csv file
                # Lets user know what gene/organism combo is being blasted next
                blastn_log.info(50 * '-')
                blastn_log.info("Here is the current list of good accession #'s: %s" % gene_list_A)
                blastn_log.info("Here is the current list of good GI #'s: %s" % gene_list_G)
                blastn_log.info(str(Organism) + " with a taxonomy id of " + str(TID) + " will be BLASTed next")
                blastn_log.info(50 * '-')

                # Create shorter variables for organism name and taxonomy id.
                Org = str(Organism)
                TAX = str(TID)

# ------------------------------------------------------------------------------
                # Standalone NCBI BLAST using local refseq_rna database on MCSR

# DEPRECATED
                # Before appending to the gene list, only append
                # the human accession # that was the input for the
                # homo sapiens list.
                # if TAX == '9606':
                #     gene_list_A.append(str(Accession[2]))
                #     BA_count = BA_count + 1
                #     cmd = "blastdbcmd -entry " + str(Accession[2]) + " -db refseq_rna -outfmt %f -out temp.fasta"
                #     getgi = "blastdbcmd -entry " + str(Accession[2]) +" -db refseq_rna -outfmt %g"
                #
                #     cmdstatus = subprocess.call([getgi], shell=True)
                #
                #     if cmdstatus == 0:  # Command was successful.
                #         pass  # Continue through script.
                #
                #     else:  # Unsuccessful. Stdout will be '1'. Entry was not found.
                #         unfound_genes.append(str(Accession[1]))  # add unfound gene to list
                #         unfound_accs.append(str(Accession[2]))
                #         unfound_tier.append(str(Accession[0]))
                #         blastn_log.error("Entry not found: %s" % str(Accession[2]))  # Log it.
                #         os.system("rm -r temp.fasta")  # Remove that temp fasta file.
                #
                #         # Remove the entry/gene/accession not found in the db from the lists.
                #         gene_list_A.remove(str(Accession[1]))
                #         gene_list_G.remove(str(Accession[1]))
                #         gene_list_TREF.remove(str(Accession[1]))
                #         gene_list_TRBF.remove(str(Accession[1]))
                #         break  # Continue to the next gene.
                #
                #     gi = subprocess.check_output([getgi], shell=True)
                #
                #     gi = gi.strip()
                #     gi = gi.decode('utf-8')
                #     gi = str(gi)
                #     gi = gi.replace("'", "")
                #     gene_list_G.append(gi)
                #     GI_count = GI_count + 1
                #
                #     blastn_log.info("----%s seconds----" % (time.time() - start_time))
                #     timer = str(time.time() - start_time)
                #     gene_list_TREF.append('')
                #     gene_list_TRBF.append(timer)
                #     break
# DEPRECATED
                else:
                    # Create a temporary fasta file since the blastn command needs a sequence file as input.
                    cmd = "blastdbcmd -entry " + str(Accession[2]) +" -db refseq_rna -outfmt %f -out temp.fasta"

                    # Call command and get the status of the command
                    status = subprocess.call([cmd], shell=True)

                    if status == 0:  # Command was successful.
                        pass  # Continue through script.

                    else:  # Unsuccessful. Stdout will be '1'. Entry was not found or some other error
                        blastn_log.error("There was an error with : %s" % str(Accession[2]))  # Log it.
                        break

                    # Create/Open a XML file that stores BLAST data for a particular Organism.
                    # By opening for writing, we can overwrite already existing xml files.
                    save_file = open("%s_%s.xml" % (Accession[1], Org), "w")

                    # Create a copy of the gi list file per taxonomy id to be used in blast
                    os.system("cp " + h + "/data/gi-lists/" + TAX + "gi " + TAX + "gi")

                    # Use Biopython's NCBIBlastnCommandline tool
                    result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna",
                                                           strand="plus", evalue=0.001,  # DONT GO LOWER
                                                           out="%s_%s.xml" % (Accession[1], Org),
                                                           outfmt=5, gilist=TAX + "gi",
                                                           max_target_seqs=10, task="blastn")
                    stdout_str, stderr_str = result_handle1()
                    #blastn_log.info(result_handle1)  # log the result handle as a check.

                    # Remove the gi list obinary file from the current directory
                    os.remove(TAX + "gi")
                    blastn_log.info(TAX + "gi file has been deleted." + "\n")

                    # Remove the temp.fasta file in the directory
                    os.remove("temp.fasta")
                    blastn_log.info("The temp.fasta file has been deleted." + "\n")
                    blastn_log.info("%s_%s.xml" % (Accession[1],Org) + " is being parsed." + "\n")

# ------------------------------------------------------------------------------
                    # Open the saved XML files above in order to sort through the BLAST Data
                    with open("%s_%s.xml" % (Accession[1], Org)) as result_handle2:
                        blast_qresult = SearchIO.read(result_handle2, 'blast-xml')
                        mapped_qresult = blast_qresult.hit_map(map_func)
                        for hit in mapped_qresult:
                            for hsp in hit.hsps:
                                blastn_log.info(hit.id1)
                                blastn_log.info(hit.id2)
                                blastn_log.info(hsp.bitscore_raw)
                                blastn_log.info(hit.description + '\n')

# ------------------------------------------------------------------------------
                                # Find the highest scoring hit for each gene
                                if hsp.bitscore_raw > maximum:
                                    # If the gene is a pseudogene then go the the next hit
                                    if "xr" in str(hit.id.lower()):
                                        blastn_log.info("Encountered a predicted(X*_00000) "
                                                        "non-coding (XR_000000)(ncRNA) RefSeq.  Moving to the next hit.")
                                        break

                                    else:
                                        # If the gene is acceptable then add it to the gene list
                                        maximum = hsp.bitscore_raw
                                        if x.lower() in hit.description.lower():
                                            Best_Accession = hit.id1
                                            Best_GI = hit.id2
                                        else:
                                            Best_Accession = hit.id1.lower()
                                            Best_GI = hit.id2
                                        blastn_log.info(Best_Accession + ' ' + Best_GI + " has the higlhest bitscore!!!!")

# ------------------------------------------------------------------------------
                                        gene_list_A.append(Best_Accession)
                                        gene_list_G.append(Best_GI)
                                        BA_count = BA_count + 1
                                        GI_count = GI_count + 1
                                        Best_Accession = ''
                                        Best_GI = ''

# ------------------------------------------------------------------------------
                        # blastn_log.info out the counts to ensure the user that everything is in order
                        blastn_log.info("BA_count = %s" % BA_count)
                        blastn_log.info("file_count = %s" % file_count)
                        # If the hit table runs out of hits then there will be no best
                        # accession number, so leave a blank in the gene_list
                        if file_count != BA_count:
                            gene_list_A.append('')
                            BA_count = file_count
                        if file_count != GI_count:
                            gene_list_G.append('')
                            GI_count = file_count

                        blastn_log.info("----%s seconds----" % (time.time() - start_time))
                        timer = str(time.time() - start_time)
                        gene_list_TREF.append('')
                        gene_list_TRBF.append(timer)

# ------------------------------------------------------------------------------
    # Add lists to csv files. If gene was not found in blastdb, add it to separate
    # file. If the gene was in the blast db, add the top hits/accessions list/
    if str(Accession[1]) in unfound_genes:  # If the gene is found in the list of unfound genes
        os.chdir(h)
        os.chdir(output)

        # Save the list of unfound genes, accessions, tiers to a csv file
        ug = pd.DataFrame(unfound_genes, dtype=str)  # create dataframe for unfound gene
        ua = pd.DataFrame(unfound_accs, dtype=str)  # create dataframe for unfound accession
        ut = pd.DataFrame(unfound_tier, dtype=str)  # create dataframe for unfound tier

        dframes = [ut, ug, ua]  # add the dataframes to a list
        norecord = pd.concat(dframes, axis=1)  # turn the dataframes into 1 dataframe

        if not os.path.isfile('genes_not_in_blastdb.csv'):  # create file if not existing
            norecord.to_csv('genes_not_in_blastdb.csv', header=None, index=False)

        else: # if csv file exists, append to it.
            norecord.to_csv('genes_not_in_blastdb.csv', mode='a', header=None, index=False)

    else:  # If the gene is present in the database, save the lists to files.
        os.chdir(h)
        os.chdir(output)
# TODO-ROB:  ERROR with extra line
        # Open the Master_Accession, Master_GI, and TREF File. Add the full gene list.
        with open('Master_Accession_File.csv', 'a', newline='') as csvfile:
            gene_row = csv.writer(csvfile, dialect='excel')
            gene_row.writerow(gene_list_A)
            blastn_log.info(check)  # Check

        with open('Master_GI_File.csv', 'a', newline='') as csvfile:
            gene_row = csv.writer(csvfile, dialect='excel')
            gene_row.writerow(gene_list_G)
            blastn_log.info(check)  # Check

        with open("Time_Record_ExistingFiles.csv", 'a', newline='') as csvfile:
            gene_row = csv.writer(csvfile, dialect='excel')
            gene_row.writerow(gene_list_TREF)
            blastn_log.info(check)  # Check

        with open("Time_Record_BLASTingFiles.csv", 'a', newline='') as csvfile:
            gene_row = csv.writer(csvfile, dialect='excel')
            gene_row.writerow(gene_list_TRBF)
            blastn_log.info(check)  # Check

# ------------------------------------------------------------------------------
# Close open files and announce the completion of the script.
end_of_script = time.time()  # End time of script/blasting

f.close()
#save_file.close()
csvfile.close()

blastn_log.info("This script took {0} minutes to complete.".format((end_of_script - start_of_script)/60))
blastn_log.info("This blastn portion of the script has completed. Check your output! ✓" + 2*("\n"))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Post Blast Analysis
# Analyze the accession file to check for missing genes

# Set up the post blastn analysis logger and log file
post_blast_log = log.getLogger('Post Blast Analysis')
# ------------------------------------------------------------------------------
# TODO-ROB:  Add post blast analysis to CompGenAnalysis
# Read output accessions file and create dictionaries for data types.
# Change to home directory
os.chdir(h)
post_blast_log.info('------------------------------------------------------------------------------')
post_blast_log.info('Post Blast Analysis has begun.' + '\n')

accession_data = OrthologAnalysis(acc_file='Master_Accession_File.csv')
missing_gene = accession_data.missing_dict['genes']
missing_orgs = accession_data.missing_dict['organisms']
orgs = accession_data.org_list

# # Read in blast output file
# original = pd.read_csv('Master_Accession_File.csv', dtype=str)
#
# # Save the master accession file to the processed folder
# final_csv = processed + 'Master_Accession_File_Final.csv'
#
# if os.path.exists(processed) == True:
#     os.system('rm -r ' + final_csv)
#     original.to_csv(processed + 'Master_Accession_File_Final.csv', dtype=str, index=False)
#     pass
# else:
#     original.to_csv(processed + 'Master_Accession_File_Final.csv', dtype=str, index=False)
#     pass
#
# # Change the index to gene column
# original = original.set_index('Gene')

# # Read in organisms file and create organisms list
# orgs = pd.read_csv(initial_data + 'organisms.csv', index_col=False, dtype=str, header=None)
# orglist = orgs.replace(to_replace=' ', value='_', regex=True)
# orglist = list(orglist[0])

# no_acc = {}  # Dictionary for no accessions by organism
# num_missing = {}  # Dictionary for # of missing genes by organism
#
# maf = original.isnull()  # Returns blanks as 'True'

# # Use a forloop to create a dictionary of missing genes by organism
# for org in orglist:
#     x = maf.ix[:, org]
#     y = org + "keylist"
#     y = []
#     for key, item in x.items():
#         if item == True:
#             y.append(key)
#             no_acc[org] = y
#             n = len(y)
#             num_missing[org] = n
#             continue

# Get number of missing organisms per gene
# miss_per_gene = maf.sum(axis=1)

# ------------------------------------------------------------------------------
# Create and write dictionaries & dataframes to excel file
if missing_gene['count'] <= 0 and missing_orgs['count'] <= 0:
    # Log that the blast had full coverage
    post_blast_log.info('There are no missing accession numbers for any gene or organism.')
    post_blast_log.info('Post blastn analysis is complete.')

else:
    post_blast_log.info('There are missing accessions. This data will be written an excel file.')

    # Set up the excel file
    excel_file = pd.ExcelWriter(processed + 'karg_missing_genes_data.xlsx')

    # This is the data frame for the names of missing genes by organism
    frame = pd.DataFrame.from_dict(no_acc, orient='index', dtype=str)
    frame = frame.transpose()
    frame.to_excel(excel_file, sheet_name="Missing Genes by Org", index=False)

    # This is the data frame for the number of missing genes by organism
    frame2 = pd.DataFrame.from_dict(num_missing, orient='index')
    frame2.to_excel(excel_file, sheet_name="# of Genes missing by Org", header=False)

    # This is the data frame for the number of missing ORGANISMS by gene
    frame3 = pd.DataFrame(miss_per_gene)
    combine = [original.Tier, frame3]  # Combine the tier column and data column
    frame3 = pd.concat(combine, axis=1) # Concatenate the dataframes
    frame3.to_excel(excel_file, sheet_name="# of Orgs missing by Gene", header=False)

    # Save the excel file
    excel_file.save()
    post_blast_log.info('Your file, karg_missing_genes_data.xlsx, has been created and saved.')

    post_blast_log.info('Post blastn analysis is complete.')
