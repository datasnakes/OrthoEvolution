import os
import sys
from Bio import SearchIO        # Used for parsing and sorting XML files.
import csv                      # Used for dealing with CSV files
from Bio.Blast.Applications import NcbiblastnCommandline  # Used for Local Blasting.
import time                     # Used to delay when dealing with NCBI server errors

# Define the map function for hit id's.
def map_func(hit):
    hit.id1 = hit.id.split('|')[3]
    hit.id2 = hit.id.split('|')[1]
    hit.id = hit.id[:-2]
    return hit

# Create a variable used to check the output
check = 1
# Create a variable used to check the processing time
start_time = time.time()

# Manage directories
y = '/home/ums/r2295/bin/Orthologs-Project/'  # Home Directory
home = y
z = '/ptmp/r2295/bin/Blast_Output/'
host = z
os.chdir(home)  # Directory Change: Home directory
os.listdir()  # Make a list of files

# Open a .csv file that contains 1 column of organism names
# Make a list of organisms and use it for column headers in the master file
org_list = []   # Initialize a list of organisms
org_list.append('')
o = open('Organisms.csv')
file2 = csv.reader(o)
for org in file2:    # Format a list of organisms
    org = str(org)
    org = org.replace("'", "")
    org = org.replace("[", "")
    org = org.replace("]", "")
    org = org.replace(" ", "_")
    org_list.append(org)
print(org_list)
# input("Ok?")

# Open a .csv file that contains 1 column of taxonomy id's
# Make a list of the tax id's and use them for local NCBI blast
# Use Taxonomy ID list to furher customize Blast search.
ID_list = []   # Initialize a list of taxidanisms
ID_list.append('')
taxid = open('taxids.csv')
file4 = csv.reader(taxid)
for ID in file4:    # Format a list of taxonomy id's
    ID = str(ID)
    ID = ID.replace("'", "")
    ID = ID.replace("[", "")
    ID = ID.replace("]", "")
    ID_list.append(ID)
print(ID_list)
# input("Ok?")

os.chdir(z) # Change or move to the 'host' directory
output_dir_list = os.listdir()  # Make a list of files

#############################################################################

# Check to see if the master ACCESSION file is in the home directory,
# and then either add the header or count the number of rows that already
# exist in order to skip to the most recently called gene.

if 'Master_Accession_File.csv' in output_dir_list:
    MAF = open('Master_Accession_File.csv')
    MAF = csv.reader(MAF)
    row_count_A = sum(1 for row in MAF) - 1
    print("row_count_A: " + str(row_count_A))
    # input("Is this an ok row_count for the Accession File?")
else:
    MAF = open('Master_Accession_File.csv', 'w', newline='')
    org_row = csv.writer(MAF, dialect='excel')
    org_row.writerow(org_list)
    MAF.close()
    row_count_A = 0
    print("2_A")

# Check to see if the master GI file is in the home directory,
# and then either add the header or count the number of rows that already
# exist in order to skip to the most recently called gene

if 'Master_GI_File.csv' in output_dir_list:
    MGF = open('Master_GI_File.csv')
    MGF = csv.reader(MGF)
    row_count_G = sum(1 for row in MGF) - 1
    print("row_count_G: " + str(row_count_G))
    # input("Is this an ok row_count for the GI File?")
else:
    MGF = open('Master_GI_File.csv', 'w', newline = '')
    org_row = csv.writer(MGF, dialect='excel')
    org_row.writerow(org_list)
    MGF.close()
    row_count_G = 0
    print("2_G")

# Check to see if the Time_Record_ExistingFiles file is in the home directory,
# and then either add the header or count the number of rows that already exist
# in order to skip to the most recently called gene

if 'Time_Record_ExistingFiles.csv' in output_dir_list:
    TREF = open('Time_Record_ExistingFiles.csv')
    TREF = csv.reader(TREF)
    row_count_TREF = sum(1 for row in TREF) - 1
    print("row_count_TREF: " + str(row_count_TREF))
    # input("Is this an ok row_count for the Time_Record_ExistingFiles File?")
else:
    TREF = open('Time_Record_ExistingFiles.csv', 'w', newline='')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(org_list)
    TREF.close()
    row_count_TREF = 0
    print("2_TREF")

# Check to see if the Time_Record_BLASTingFiles file is in the home directory,
# and then either add the header or count the number of rows that already exist
# in order to skip to the most recently called gene

if 'Time_Record_BLASTingFiles.csv' in output_dir_list:
    TRBF = open('Time_Record_BLASTingFiles.csv')
    TRBF = csv.reader(TRBF)
    row_count_TRBF = sum(1 for row in TRBF) - 1
    print("row_count_TRBF: " + str(row_count_TRBF))
    # input("Is this an ok row_count for the Time_Record_BLASTingFiles File?")
else:
    TREF = open('Time_Record_BLASTingFiles.csv', 'w', newline='')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(org_list)
    TREF.close()
    row_count_TREF = 0
    print("2_TRBF")

# Output Check
    # input("Is this an ok row_count for the Accession File, GI File, TREF file, and TRBF file?")

# For writing Accessions/GIs to the file
gene_list_A = []
gene_list_B = []
file_count = 0
BA_count = 0
GI_count = 0

#############################################################################

os.chdir(home)

# blast the Homo Sapiens Accession number & filter results via the Organism
# The 1st 'for' loop parses through the individual genes in our data files

f = open('Homo_sapiens_Accession.csv')  # 1st column - gene names;  2nd column - Human Acc. No.s
file1 = csv.reader(f)
m = open('Macaca_mulatta_Accession.csv')  # 1st column - gene names;  2nd column - Rhesus Acc. No.s
file3 = csv.reader(m)
os.chdir(z)

Acc_count = 0  # The Accession count. Lists start at 0.

for Accession in file1:
    Acc_count = Acc_count + 1

    for monkey in file3:
        Rhesus = str(monkey[1])  # Rhesus Accession number
        break

# Skip over genes that have already been recorded in the Master_File
# Comment this part out if you need to rebuild any of the master files
#       if Acc_count <= row_count_A: # something was deleted here -> (and....)
#       input("This gene %s has been accounted for in the Master_File" % Accession[0])
#       continue

#  Begin listing I/O information
#  ***********************************************************************************************************
#  ***********************************************************************************************************

    print(50 * '*')
    print(50 * '*')
    print('The following contains information about the blast input and output:')
    print('Human Accesion: %s' % Accession[1])
    print('Rhesus Accession: %s' % Rhesus)
    print("Gene of Interest: %s" % Accession[0])

    os.chdir(host)
    print("Current Directory: " + str(os.getcwd()))  # Check the Python Shell for proper directory

    gene_list_A = []
    gene_list_G = []
    gene_list_TREF = []
    gene_list_TRBF = []

    gene_list_A.append(str(Accession[0]))  # Gene name for rows
    #gene_list_A.append(str(Accession[1]))  # Human Accession
    #gene_list_A.append(str(Rhesus))  # Rhesus Accession

    gene_list_G.append(str(Accession[0]))
    #gene_list_G.append(str(Accession[1]))
    #gene_list_G.append(str(Rhesus))

    gene_list_TREF.append(str(Accession[0]))
    #gene_list_TREF.append(str(Accession[1]))
    #gene_list_TREF.append(str(Rhesus))

    gene_list_TRBF.append(str(Accession[0]))
    #gene_list_TRBF.append(str(Accession[1]))
    #gene_list_TRBF.append(str(Rhesus))


# Create a folder (for target genes) unless it already exists.
# If it does exist, then change to that folder to see
# what XML files (for target gene/organisms) are present.

    x = str(Accession[0])  # Used for making/changing directories

    try:
        os.mkdir('%s' % (x))
        print("Directory Created: %s" % (x))
        os.chdir(host + x)
        print("Current Directory: " + str(os.getcwd()))
        print("\n")

    except FileExistsError:
        print("Directory already exists: %s" % (x))
        os.chdir(host + x)
        print("Current Directory: " + str(os.getcwd()))

    print(os.getcwd())  # Check the Python Shell for proper output
    print(50 * '*')
    print((50 * '*') + '\n' + '\n')

#  End listing I/O information
#  ***********************************************************************************************************
#  ***********************************************************************************************************


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################


# The 2nd 'for' loop is for parsing through the individual organisms(which correspond
# to taxonomy id's) in our data files, so that very specific blast data can be
# retrieved and stored in the above directory.
    for Organism, TID in zip(org_list, ID_list):
        file_count = file_count + 1
        maximum = 0
        os.chdir(z + x) # Change the directory

# Skip the first 3 items in the org_list and update the Best Accession Count and Best GI Count
        if Organism == '':  # or Organism == 'Homo sapiens' or  Organism == 'Macaca mulatta':
            BA_count = file_count
            GI_count = file_count
            continue
        xmlfile_list = os.listdir()
        s = "%s_%s.xml" % (Accession[0], Organism)

# If a gene/organism XML file has already been established, parse through it.
        if s in xmlfile_list:
            print(Organism)
            print(org_list)
            with open("%s_%s.xml" % (Accession[0], Organism)) as result_handle2:
                blast_qresult = SearchIO.read(result_handle2, 'blast-xml')
                mapped_qresult = blast_qresult.hit_map(map_func)
                for hit in mapped_qresult:
                    for hsp in hit.hsps:
                        print(hit.id1)
                        print(hit.id2)
                        print(hsp.bitscore_raw)
                        print(hit.description + '\n')
                        if hsp.bitscore_raw > maximum:
                            if "xr" in str(hit.id.lower()):
                                print("Encountered a pseudogene.  Moving to the next hit.")
                                break
                            else:
                                maximum = hsp.bitscore_raw
                                if x.lower() in hit.description.lower():
                                    Best_Accession = hit.id1
                                    Best_GI = hit.id2
                                else:
                                    Best_Accession = hit.id1.lower()
                                    Best_GI = hit.id2
                                print(Best_Accession + ' ' + Best_GI + " has the higlhest bitscore!!!!")
                                print(hsp.hit)
                                gene_list_A.append(Best_Accession)
                                gene_list_G.append(Best_GI)
                                BA_count = BA_count + 1
                                GI_count = GI_count + 1
                print("BA_count = %s" % BA_count)
                print("GI_count = %s" % GI_count)
                print("file_count = %s" % file_count)
                if file_count != BA_count:
                    gene_list_A.append('')
                    BA_count = file_count
                if file_count != GI_count:
                    gene_list_G.append('')
                    GI_count = file_count
                print("----%s seconds----" % (time.time() - start_time))
                timer = str(time.time() - start_time)
                gene_list_TREF.append(timer)
                gene_list_TRBF.append('')

            print(s + " already exists.  Next organism.")
            # file_count = file_count + 1
            continue

# If the gene/organism xml file hasn't already been established,
# then blast the Homo sapiens accession number,
# and store the blast report in an XML file.
        else:
            # Lets user know what gene/organism combo is being blasted next
            print(50 * '$')
            print("Here is the current list of good accession #'s: %s" % gene_list_A)
            print("Here is the current list of good GI #'s: %s" % gene_list_G)
            print(str(Organism) + " with a taxonomy id of " + str(TID) + " will be BLASTed next")
            print(50 * '$')

            Org = str(Organism)
            TAX = str(TID)

            # Create/Open a XML file that stores blast data for a particular Organism.
            # By opening for writing, we can overwrite already existing xml files.

            save_file = open("%s_%s.xml" % (Accession[0],Org), "w")

        ### Standalone NCBI blast using local refseq_rna database on MCSR ###

            # Create a copy of the gi list file per taxonomy id to be used in blast
            os.system("cp /work2/vallender/Projects/GPCR-Orthologs/DocsAndFiles/gi_lists/" + TAX + "gi " + TAX + "gi")

            # Create a temporary fasta file since the blastn command needs a sequence file as input.
            os.system("blastdbcmd -entry "+ str(Accession[1]) +" -db refseq_rna -outfmt %f -out temp.fasta")
            print(open('temp.fasta', 'r').read())

            # Use Biopython's NCBIBlastnCommandline tool
            result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna", strand="plus", evalue=0.0005,
                                                   out="%s_%s.xml" % (Accession[0], Org), word_size=20, best_hit_score_edge=0.05,
                                                   outfmt=5, gilist=TAX + "gi", max_target_seqs=10, task="blastn")
            stdout_str, stderr_str = result_handle1()
            print(result_handle1)  # Print the result handle as a check.

            # Remove the gi list obinary file from the current directory
            os.remove(TAX + "gi")
            print("\n" + TAX + "gi file has been deleted.")
            time.sleep(.5)

            # Remove the temp.fasta file in the directory
            os.remove("temp.fasta")
            print("\n" "The temp.fasta file has been deleted.")
            time.sleep(.5)
            print("\n" + "%s_%s.xml" % (Accession[0],Org) + " is being parsed." + "\n")
            time.sleep(.2)

# Open the saved XML files above in order to sort through the blast Data
            with open("%s_%s.xml" % (Accession[0], Org)) as result_handle2:
                blast_qresult = SearchIO.read(result_handle2, 'blast-xml')
                mapped_qresult = blast_qresult.hit_map(map_func)
                for hit in mapped_qresult:
                    for hsp in hit.hsps:
                        print(hit.id1)
                        print(hit.id2)
                        print(hsp.bitscore_raw)
                        print(hit.description + '\n')

# Find the highest scoring hit for each gene
                        if hsp.bitscore_raw > maximum:
                            # If the gene is a pseudogene then go the the next hit
                            if "xr" in str(hit.id.lower()):
                                print("Encountered a pseudogene.  Moving to the next hit.")
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
                                print(Best_Accession + ' ' + Best_GI + " has the higlhest bitscore!!!!")
                                gene_list_A.append(Best_Accession)
                                gene_list_G.append(Best_GI)
                                BA_count = BA_count + 1
                                GI_count = GI_count + 1
                                Best_Accession = ''
                                Best_GI = ''

                # Print out the counts to ensure the user that everything is in order
                print("BA_count = %s" % BA_count)
                print("file_count = %s" % file_count)
                # If the hit table runs out of hits then there will be no best accession number, so leave a blank in the gene_list
                if file_count != BA_count:
                    gene_list_A.append('')
                    BA_count = file_count
                if file_count != GI_count:
                    gene_list_G.append('')
                    GI_count = file_count

                print("----%s seconds----" % (time.time() - start_time))
                timer = str(time.time() - start_time)
                gene_list_TREF.append('')
                gene_list_TRBF.append(timer)

    os.chdir(z)

# Open the Master_Accession, Master_GI, and TREF File. Add the full gene list.
    with open('Master_Accession_File.csv', 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_A)

        print(check)  # Check

    with open('Master_GI_File.csv', 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_G)

        print(check)  # Check

    with open("Time_Record_ExistingFiles.csv", 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_TREF)

        print(check)  # Check

    with open("Time_Record_BLASTingFiles.csv", 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_TRBF)

        print(check)  # Check

    o.close()  # File Handling

    os.chdir(y)  # Directory Change: Home directory

# Close open files and announce the completion of the script.
f.close()
m.close()
o.close()
save_file.close()
print("\n")
sys.exit("This script has completed. Check your output! ✓")