#################################################################################################

#title           :HTR1A_1.0.py
#description     :This is PART A (Generating Accessions) of the Orthologs Project.  See Outline.
#       
#author          :Rob Gilmore RAG
#date            :2/3/2016
#version         :0.1
#notes           :
#python_version  :3.4.3; BioPython

##### IT IS IMPORTANT TO NOTE THAT THIS CODE IS DESIGNED TO USE NCBI's SERVERS IN ORDER TO blast
##### THE DESIRED SEQUENCES.  WHEN ACTUALLY IMPLEMENTING THIS CODE, WE WILL BE USING THE
##### SUPERCOMPUTER, WHICH MAY OR MAY NOT CONTAIN A LOCAL VERSION OF blast.  IF IT THE
##### SUPERCOMPUTER HAS A LOCAL VERSION OF blast THEN THIS CODE WILL HAVE TO CHANGE.

#################################################################################################

from Bio import SeqIO           ## Will be used for alignment (PART B)
from Bio import AlignIO         ## Will be used for alignment (PART B)
from Bio import SearchIO        ## Used for parsing and sorting XML files.
from Bio.Blast import NCBIWWW   ## Used for BLASTING over the internet
from Bio.Blast import NCBIXML   ## Used for parsing XML files.
from Bio.Blast.Applications import NcbiblastnCommandline  ##  Used for Local Blasting.
import os                       ## Used for organizing directories
import csv                      ## Used for dealing with CSV files
import time                     ## Used to delay when dealing with NCBI server errors
import xlsxwriter               ## Used to format output.

def map_func(hit):
    hit.id1 = hit.id.split('|')[3]
    hit.id2 = hit.id.split('|')[1]
    hit.id = hit.id[:-2]
    return hit

## Create a variable used to check the output
check = 1
## Create a variable used to check the processing time
start_time = time.time()



y = 'X:\\Bioinformatics\\Ortholog Project\\blast Input\\' ## Home Directory
home = y
z = y + '\\blast Output\\'
host = z
os.chdir(home) ## Directory Change: Home directory
os.listdir() ## Make a list of files


## Open a CSV file that contains 1 column of organism names
## Make a list of organisms and use it for column headers in the master file
org_list = []  ## Initialize list of Organisms
org_list.append('')
o = open('Organisms.csv')
file2 = csv.reader(o)
for org in file2:  ##  Format a list of organisms
    org = str(org)
    org = org.replace("'", "")
    org = org.replace("[", "")
    org = org.replace("]", "")
    org_list.append(org)
print(org_list)
os.chdir(z)
output_dir_list = os.listdir() ## Make a list of files
##########################################################################################################


##  Check to see if the master ACCESSION file is in the home directory, and then either add the header or count
##  the number of rows that already exist in order to skip to the most recently called gene
if 'Master_Accession_File.csv' in output_dir_list:
    MAF = open('Master_Accession_File.csv')
    MAF = csv.reader(MAF)
    row_count_A = sum(1 for row in MAF) - 1
    print("row_count_A: " + str(row_count_A))
    input("Is this an ok row_count for the Accession File?")
else:
    MAF = open('Master_Accession_File.csv', 'w', newline = '')
    org_row = csv.writer(MAF, dialect='excel')
    org_row.writerow(org_list)
    MAF.close()
    row_count_A = 0
    print("2_A")

    
##  Check to see if the master GI file is in the home directory, and then either add the header or count
##  the number of rows that already exist in order to skip to the most recently called gene
if 'Master_GI_File.csv' in output_dir_list:
    MGF = open('Master_GI_File.csv')
    MGF = csv.reader(MGF)
    row_count_G = sum(1 for row in MGF) - 1
    print("row_count_G: " + str(row_count_G))
    input("Is this an ok row_count for the GI File?")
else:
    MGF = open('Master_GI_File.csv', 'w', newline = '')
    org_row = csv.writer(MGF, dialect='excel')
    org_row.writerow(org_list)
    MGF.close()
    row_count_G = 0
    print("2_G")

##  Check to see if the Time_Record_ExistingFiles file is in the home directory, and then either add the header or count
##  the number of rows that already exist in order to skip to the most recently called gene
if 'Time_Record_ExistingFiles.csv' in output_dir_list:
    TREF = open('Time_Record_ExistingFiles.csv')
    TREF = csv.reader(TREF)
    row_count_TREF = sum(1 for row in TREF) - 1
    print("row_count_TREF: " + str(row_count_TREF))
    input("Is this an ok row_count for the Time_Record_ExistingFiles File?")
else:
    TREF = open('Time_Record_ExistingFiles.csv', 'w', newline = '')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(org_list)
    TREF.close()
    row_count_TREF = 0
    print("2_TREF")

##  Check to see if the Time_Record_BLASTingFiles file is in the home directory, and then either add the header or count
##  the number of rows that already exist in order to skip to the most recently called gene
if 'Time_Record_BLASTingFiles.csv' in output_dir_list:
    TRBF = open('Time_Record_BLASTingFiles.csv')
    TRBF = csv.reader(TRBF)
    row_count_TRBF = sum(1 for row in TRBF) - 1
    print("row_count_TRBF: " + str(row_count_TRBF))
    input("Is this an ok row_count for the Time_Record_BLASTingFiles File?")
else:
    TREF = open('Time_Record_BLASTingFiles.csv', 'w', newline = '')
    org_row = csv.writer(TREF, dialect='excel')
    org_row.writerow(org_list)
    TREF.close()
    row_count_TREF = 0
    print("2_TRBF")
    
## Output Check
input("Is this an ok row_count for the Accession File, GI File, TREF file, and TRBF file?")

##  For writing Accessions/GI to the file

gene_list_A = []
gene_list_B = []
file_count = 0
BA_count = 0
GI_count = 0


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################

os.chdir(home)
## blast the Homo Sapien Accession number and filter the results via the Organism list
## The 1st 'for' loop is for parsing throught the individual genes in our data files
f = open('Homo_sapien_Accession.csv')  ## 1st column - gene names;  2nd column - Human Acc. No.s
file1 = csv.reader(f)
m = open('Macaca_mulatta_Accession.csv')  ## 1st column - gene names;  2nd column - Rhesus Acc. No.s
file3 = csv.reader(m)
os.chdir(z)

Acc_count = 0

for Accession in file1:
    Acc_count = Acc_count + 1

    for monkey in file3:
        Rhesus = str(monkey[1]) ##  Rhesus Acc. No.
        break
    
##  Skip over genes that have already been recorded in the Master_File
##  Comment this part out if you need to rebuild any of the master files
##    if Acc_count <= row_count_A: # something was deleted here -> (and....)
##        input("This gene %s has been accounted for in the Master_File" % Accession[0])
##        continue

    
##  Begin listing I/O information
##  ***********************************************************************************************************
##  ***********************************************************************************************************
    
    print('****************************************************************************************************')
    print('****************************************************************************************************')
    print('The following contains information about the blast input and output:')
    print('Human Accesion: %s' % Accession[1])
    print('Rhesus Accession: %s' % Rhesus)
    print("Gene of Interest: %s" % Accession[0])     

    os.chdir(host)
    print("Current Directory: " + str(os.getcwd())) ## Check the Python Shell for proper directory
    
    gene_list_A = []
    gene_list_G = []
    gene_list_TREF = []
    gene_list_TRBF = []
    
    gene_list_A.append(str(Accession[0]))  ## Gene name for rows
    #gene_list_A.append(str(Accession[1]))  ## Human Accession
    #gene_list_A.append(str(Rhesus))  ## Rhesus Accession
    
    gene_list_G.append(str(Accession[0]))
    #gene_list_G.append(str(Accession[1]))
    #gene_list_G.append(str(Rhesus))
    
    gene_list_TREF.append(str(Accession[0]))
    #gene_list_TREF.append(str(Accession[1]))
    #gene_list_TREF.append(str(Rhesus))

    gene_list_TRBF.append(str(Accession[0]))
    #gene_list_TRBF.append(str(Accession[1]))
    #gene_list_TRBF.append(str(Rhesus))
    
    
## Try to create a folder (for target genes) unless it already exists.  If it does then change to that folder to see
## what XML files (for target gene/organisms) are present.

    x = str(Accession[0])  ## Used for making/changing directories
    
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
        
    
    
    print(os.getcwd()) ## Check the Python Shell for proper output
    print('****************************************************************************************************')
    print('****************************************************************************************************' + '\n' + '\n')
    
##  End listing I/O information
##  ***********************************************************************************************************
##  ***********************************************************************************************************
            

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################

    
## The 2nd 'for' loop is for parsing through the individual organisms in our data files,
## so that very specific blast data can be retrieved and stored in the above directory.
    for Organism in org_list:
        file_count = file_count + 1
        maximum = 0
        
        os.chdir(z + x)
        
## Skip the first 3 items in the org_list and update the Best Accession Count and Best GI Count            
        if Organism == '': #or Organism == 'Homo sapiens' or  Organism == 'Macaca mulatta':
            BA_count = file_count
            GI_count = file_count
            continue
    
        xmlfile_list = os.listdir()
        s = "%s_%s.xml" % (Accession[0], Organism)

##  If a gene/organism XML file has already been established, then parse through it.        
        if s in xmlfile_list:
            with open("%s_%s.xml" % (Accession[0],Organism)) as result_handle2:
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

        
##  If the gene/organism xml file hasn't already been established, then blast the HS accession number,
##  and store the blast report in an XML file.
        else:  
##  Let the user know what gene/organism combination is being blasted next            
            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            print("Here is the current list of good accession #'s: %s" % gene_list_A)
            print("Here is the current list of good GI #'s: %s" % gene_list_G)
            print(str(Organism) + " will be BLASTed next")
            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            
            Org = str(Organism)
            
##  Create/Open a XML file that stores blast data for a particular Organism.
##  By opening for writing we can overwrite already existing xml files.
            
            save_file = open("%s_%s.xml" % (Accession[0],Org), "w")


##            result_handle1 = NcbiblastnCommandline(remote ="refseq_rna", out = "%s_%s.xml" % (Accession[0],Org), query = Accession[1],db = "refseq_rna", strand = "plus", evalue = 0.001, outfmt = 5, entrez_query = str(Org) +'[Organism]',)
##            print(result_handle1)
##            stdout, stderr = result_handle1()
##            input("Verdict?")
##########SET UP THE ENVIRONMENT VARIABLE TO GE THIS CODE WORKING FOR THE LOCAL blast#############

            
            
#_____________________________________________NCBI BLAST_________________________________________________________________            
            NCBI_flag = 0
            while NCBI_flag == 0:
                try:
##  1 handle, 10 hit matches, 1 target organism, 1 human orthologous gene (Accession[1])
                    result_handle1 = NCBIWWW.qblast("blastn","refseq_rna",Accession[1],entrez_query = str(Org) +'[Organism]', hitlist_size=10)
                    NCBI_flag = 1 
                    
                except:
##  NCBI servrs stop responding on occasion
                    print("NCBI servers may have timed out.  Attempting to reconnect....")
                    NCBI_flag = 0
                    time.sleep(10)
#_________________________________________________________________________________________________________________________
                    
            save_file.write(result_handle1.read()) ## Handle can only be used once after the initial blast
            save_file.close()
            result_handle1.close()

## Open the saved XML files above in order to sort through the blast Data
            with open("%s_%s.xml" % (Accession[0],Org)) as result_handle2:
                blast_qresult = SearchIO.read(result_handle2, 'blast-xml')
                mapped_qresult = blast_qresult.hit_map(map_func)
                for hit in mapped_qresult:
                    for hsp in hit.hsps:
                        print(hit.id1)
                        print(hit.id2)
                        print(hsp.bitscore_raw)
                        print(hit.description + '\n')
                        
##  Find the highest scoreing hit for each gene 
                        
                        if hsp.bitscore_raw > maximum:
##  If the gene is a pseudogene then go the the next hit
                            if "xr" in str(hit.id.lower()):
                                print("Encountered a pseudogene.  Moving to the next hit.")
                                break
                        
                            else:
##  If the gene is acceptable then add it to the gene list
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
                            
##  Print out the counts to ensure the user that everything is in order
                print("BA_count = %s" % BA_count)
                print("file_count = %s" % file_count)
##  If the hit table runs out of hits then there will be no best accession number, so leave a blank in the gene_list
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

##  Open the Master_Accession, Master_GI, and TREF File and add the full gene list.
    with open('Master_Accession_File.csv','a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_A)

        print(check) #Check

    with open('Master_GI_File.csv','a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_G)

        print(check) #Check

    with open("Time_Record_ExistingFiles.csv", 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_TREF)

        print(check) #Check
        
    with open("Time_Record_BLASTingFiles.csv", 'a', newline='') as csvfile:
        gene_row = csv.writer(csvfile, dialect='excel')
        gene_row.writerow(gene_list_TRBF)

        print(check) #Check

        
    o.close() #File Handling    

    os.chdir(y) #Directory Change: Home directory
    
        
f.close()
o.close()
save_file.close()
