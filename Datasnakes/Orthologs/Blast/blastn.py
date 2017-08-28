"""Optimized for use with local/standalone NCBI BLAST 2.6.0."""
import csv
import os
import shutil
import subprocess
import time  # Used to delay when dealing with NCBI server errors
from datetime import datetime as d
from pathlib import Path
import pandas as pd
import pkg_resources
from Datasnakes.Manager import index
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline
from Datasnakes.Orthologs.CompGenetics.ncbi_blast import BLASTAnalysis as BT
# TODO-ROB: Find packages for script timing and analysis


class BLASTn(BT):
    """Use BLASTn to search nucleotide databases using a nucleotide query.
    This class currently only works with the standalone blast.
    """

    def __init__(self, project, template=None, save_data=True, **kwargs):
        """Inherit from the BLASTing Template."""
        super().__init__(project=project, template=template, save_data=save_data, **kwargs)
        # # TODO-ROB Add taxon parameter
        # Manage Directories
        self.__home = Path(os.getcwd())
        self.__output_path = self.raw_data / Path('BLAST')  # Output directory
        self.__gi_list_path = self.__output_path / Path('gi_lists')
        self.__xml_path = self.__output_path / Path('xml')
        Path.mkdir(self.__output_path, parents=True, exist_ok=True)
        Path.mkdir(
            self.__gi_list_path /
            Path('data'),
            parents=True,
            exist_ok=True)
        Path.mkdir(self.__xml_path, parents=True, exist_ok=True)
        # # Initialize Logging
        # self.__blastn_log = LogIt.blastn()
        #df = LogIt()
        self.date_format = '%a %b %d at %I:%M:%S %p %Y'
        # self.get_time = time.time  # To get the time use 'get_time()'
        # TODO-ROB:  Add a query organism variable
        self.query_gi_dict = {}
        self.removed_genes = []
        # TODO-ROB:  Set up blast config logger, blasting logger, and post blast analysis logger
        self.blastn_log.info("These are the organisms: " + str(self.org_list))
        self.blastn_log.info("These are the genes: " + str(self.gene_list))
        self.blastn_log.info(
            "These are the taxonomy ids: \n\n\n" + str(self.taxon_ids))
        # ---------------------------------------------------------------------
        # For completed blast files
        self.complete_file = self.project + '_MAF.csv'
        self.complete_file_path = self.data / Path(self.complete_file)
        self.complete_time_file = self.project + '_TIME.csv'
        self.complete_time_file_path = self.data / Path(self.complete_time_file)

    @staticmethod
    def map_func(hit):
        """Use the map function for formatting hit id's.
        This will be used later in the script.
        """
        hit.id1 = hit.id.split('|')[3]
        hit.id2 = hit.id.split('|')[1]
        hit.id = hit.id[:-2]
        return hit

    def blast_config(self, query_align, query_organism, auto_start=False):
        """Configure everything for a BLAST.
        First the accession file, and gene list is configured.
        """
        # os.chdir(str(self.__output_path))
        self.blastn_log.info(
            '***********************************BLAST CONFIG START************ \
            ***********************\n\n\n')
        self.blastn_log.info('Configuring the accession file...')

        # Update the gene_list based on the existence of a incomplete blast
        # file
        gene_list = self.blast_file_config(self.building_file_path)
        if gene_list is not None:
            # Number of genes already BLASTed
            start = len(self.blast_human) - len(gene_list)
            # Reconfigure query_align to reflect the existing accession info
            query_align = self.blast_human[start:]
            # Reconfigure the gene_list to reflect the existing accession info
            new_gene_list = gene_list
        else:
            new_gene_list = self.gene_list

        # Create GI lists
        self.blastn_log.info(
            "Configuring GI list using the taxonomy id and the blastdbcmd tool.")
        self.gi_list_config()
        # Get GI (stdout) and query sequence (FASTA format)
        self.blastn_log.info("Generating directories.")
        self.blastn_log.info("Extracting query gi number to stdout and "
                             "query refseq sequence to a temp.fasta file from BLAST database.")
        # Iterate the query accessions numbers
        for query in query_align:
            # os.chdir(str(self.__output_path))
            gene = self.acc_dict[query][0][0]
            gene_path = self.__xml_path / Path(gene)
            org = self.acc_dict[query][0][1]
            # Create the directories for each gene
            try:
                Path.mkdir(gene_path)
                self.blastn_log.info("Directory Created: %s" % gene)
                self.blastn_log.info("\n")
                # os.chdir(gene)
            except FileExistsError:
                self.blastn_log.info("Directory already exists: %s" % gene)
                # os.chdir(gene)

            # Save sequence data in FASTA file format and print the gi number to stdout with a custom BLAST extraction
            # https://www.ncbi.nlm.nih.gov/books/NBK279689/#_cookbook_Custom_data_extraction_and_form_
            # TODO-SDH Combine these BLAST extractions???
            fmt = {
                'query': query,
                'temp fasta': str(
                    gene_path /
                    Path('temp.fasta'))}
            fasta_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %f -out {temp fasta}".format(
                **fmt)
            fasta_status = subprocess.call([fasta_setup], shell=True)
            gi_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %g".format(
                **fmt)
            gi_status = subprocess.call([gi_setup], shell=True)
            # TODO-ROB:  Add function to add the gi numbers to the dataframe/csv-file,
            # TODO-ROB: and add a check function to see if thats already there
            # Check the status of the custom blast data extraction
            if gi_status == 0 or fasta_status == 0:  # Command was successful.
                if gi_status != 0:
                    # Log it.
                    self.blastn_log.error(
                        "GI number for %s not found in the BLAST extraction" %
                        query)
                    # TODO-SDH Is this the correct move below???
                    self.blastn_log.error(
                        "Removing %s from the BLAST list..." % gene)
                    self.gene_list.remove(gene)
                    self.removed_genes.append(gene)
                    continue
                if fasta_status != 0:
                    self.blastn_log.error(
                        "FASTA sequence for %s not found in the BLAST extraction" %
                        query)
                    self.blastn_log.error(
                        "Removing %s from the BLAST list..." % gene)
                    self.gene_list.remove(gene)
                    self.removed_genes.append(gene)
                    continue
                pass
            else:
                self.blastn_log.error(
                    "FASTA sequence and GI number for %s not found in the custom BLAST extraction." %
                    query)
                self.blastn_log.error(
                    "Removing %s from the BLAST list..." %
                    gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)
                continue

            # Get the gi number from stdout, format it, and add it to the gi
            # dictionary
            gi = subprocess.check_output([gi_setup], shell=True)
            gi = gi.strip()
            gi = gi.decode('utf-8')
            gi = str(gi)
            gi = gi.replace("'", "")
            self.query_gi_dict[gene] = gi

        new_query_align = query_align
        self.blastn_log.info(
            'Configured query accession list: %s' %
            new_query_align)
        self.blastn_log.info('Configured gene list: %s\n\n\n' % new_gene_list)
        self.blastn_log.info(
            '************************************BLAST CONFIG END************************************\n\n\n')
        if auto_start is True:
            # Automatically begin BLASTING after the configuration
            self.blasting(
                genes=new_gene_list,
                query_organism=query_organism,
                pre_configured=auto_start)
        else:
            # Manually begin BLASTING and return the new gene and new query
            # lists
            return new_gene_list

    def gi_list_config(self):
        # TODO-ROB THis is for development / testing
        # TODO-ROB Add the ability to do two seperate gi configs
        """Create a gi list based on the refseq_rna database for each taxonomy id on the MCSR.
        It will also convert the gi list into a binary file which is more
        efficient to use with NCBI's Standalone Blast tools.
        """
        print('gi_list_config')
        # Directory and file handling
        cd = os.getcwd()
        os.chdir(str(self.__gi_list_path))
        taxids = self.taxon_ids
        Path.mkdir(self.__gi_list_path / Path('data'), parents=True, exist_ok=True)
        pd.Series(taxids).to_csv('taxids.csv', index=False)
        # PBS job submission using the templates
        pbs_script = 'get_gi_lists.sh'
        pbs_script_path = self.__gi_list_path / Path(pbs_script)
        py_script = 'get_gi_lists.py'
        shutil.copy(pkg_resources.resource_filename(index.__name__, pbs_script), self.raw_data)
        shutil.copy(pkg_resources.resource_filename(index.__name__, py_script), self.raw_data)
        gi_config = subprocess.check_output('qsub %s' % str(pbs_script_path), shell=True)
        gi_config = gi_config.decode('utf-8')
        print('The GI list configuration\'s JobID is %s' % gi_config)
        job_id = gi_config.replace('.sequoia', '')
        time.sleep(20)  # Wait for the job to be queued properly
        while True:
            out, err = subprocess.Popen('qsig -s SIGNULL %s' % job_id, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            print("Waiting...")
            time.sleep(30)
            if err.decode('utf-8') == 'qsig: Request invalid for state of job %s.sequoia\n' % job_id:
                print('The blast config is in MCSR\'s queue.  Waiting...')
                continue
            else:
                print('out:', out)
                print('err:', err)
                break
        os.chdir(cd)

    def blast_file_config(self, file):
        """Create or use a blast configuration file.
        This function configures different files for new BLASTS.
        It also helps recognize whether or not a BLAST was terminated
        in the middle of the dataset.  This removes the last line of
        the accession file if it is incomplete.
        """
        global ending
        output_dir_list = os.listdir(
            self.__output_path)  # Make a list of files
        # If the file exists then make a gene list that picks up from the last
        # BLAST
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

                self.blastn_log.info(
                    "The last row looks like: \n\t%s\n\t%s\n" %
                    (self.header, ending))
                self.blastn_log.info(
                    "The BLAST ended on the following query: \n%s" %
                    ncbi)
                if len(ending) < len(self.header):
                    self.blastn_log.info(
                        "Restarting the BLAST for the previous gene...")
                    count = count - 1
                # The continued gene list starts with the previous gene.
                continued_gene_list = list(
                    x for i, x in enumerate(
                        self.gene_list, 1) if i > count)
            return continued_gene_list
        # If the file doesn't exist return nothing
        else:
            self.blastn_log.info("A new BLAST started at %s" % self.get_time())
            return None

    def blast_xml_parse(self, xml_file, gene, organism):
        """Parse the XML file created by the BLAST."""
        global gi, raw_bitscore
        self.blastn_log.info(
            "Parsing %s to find the best accession number." %
            xml_file)
        maximum = 0
        file_path = str(Path(self.__xml_path) / Path(gene) / Path(xml_file))
        with open(file_path, 'r') as blast_xml:
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
            self.blastn_log.info(
                "The best accession has been selected from the BLAST xml record.")
            self.blastn_log.info("Accession:  %s" % accession)
            self.blastn_log.info(f"GI number: {gi}")
            self.blastn_log.info("Raw bitscore: %s" % raw_bitscore)
            self.blastn_log.info("Description: %s" % description)
            self.add_accession(gene, organism, accession)

    def blasting(self, genes=None, query_organism=None, pre_configured=False):
        """Configure the BLAST."""
        if pre_configured is False:
            query = self.df[query_organism].tolist()
            genes = self.blast_config(
                query_align=query,
                query_organism=query_organism,
                auto_start=True)
        elif pre_configured is True:
            genes = genes

        self.blastn_log.info(
            "------------------------------------------------------------------")
        self.blastn_log.info(
            f"The script name is {os.path.basename(__file__)}")
        self.blastn_log.info(
            'The script began on {}'.format(str(
                d.now().strftime(
                    self.date_format))))
        self.blastn_log.info(
            "------------------------------------------------------------------")

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
                # TODO-ROB change the __gi_list_path to current path + 'data'
                taxon_id = self.taxon_dict[organism]
                taxon_gi_file = str(taxon_id) + "gi"
                taxon_gi_path = self.__gi_list_path / \
                    Path('data') / Path(taxon_gi_file)
                # taxgi_dest_path = gene_path / Path(taxon_gi_file)
                if xml in files:
                    self.blast_xml_parse(xml_path, gene, organism)
                else:
                    self.blastn_log.warning(
                        "\n\n\n*******************BLAST START*******************")
                    start_time = self.get_time()
                    self.blastn_log.info("The start time is %s" % start_time)
                    self.blastn_log.info(
                        "The current gene is %s (%s)." %
                        (gene, self.tier_dict[gene]))
                    self.blastn_log.info(
                        "The current organisms is %s (%s)." %
                        (organism, taxon_id))

                    with open(xml_path, 'w') as blast_xml:
                        # TODO-ROB: For multiporcessing copy gi lists, but for regular processing just use the one.
                        # Create a copy of the gi list file per taxonomy id to be used in blast
                        # fmt = {'src': str(taxon_gi_path), 'dst': str(taxgi_dest_path)}
                        # os.system("cp {src} {dst}".format(**fmt))

                        # Set up blast parameters
                        taxon_gi_path = str(taxon_gi_path)
                        query_seq_path = str(gene_path / Path('temp.fasta'))

                        # Use Biopython's NCBIBlastnCommandline tool
                        result_handle1 = NcbiblastnCommandline(query=query_seq_path, db="refseq_rna",
                                                               strand="plus", evalue=0.001,  # DONT GO LOWER
                                                               outfmt=5, gilist=taxon_gi_path,
                                                               max_target_seqs=10, task="blastn")
                        # Capture blast data
                        stdout_str, stderr_str = result_handle1()
                        blast_xml.write(stdout_str)
                        end_time = self.get_time()
                        elapsed_time = end_time - start_time
                        self.blastn_log.info("%s was create." % blast_xml.name)
                        self.blastn_log.info("The end time is %s." % end_time)
                        self.blastn_log.info(
                            "The BLAST took %s." %
                            elapsed_time)
                    self.blastn_log.warning(
                        "********************BLAST END********************\n\n\n")
                    self.add_blast_time(gene, organism, start_time, end_time)
                    self.blast_xml_parse(xml, gene, organism)

                # If the BLAST has gone through all orthologs then create a
                # master accession file.
                if gene == genes[-1] and organism == self.org_list[-1]:
                    shutil.copy(str(self.building_file_path),
                                str(self.complete_file_path))
                    shutil.copy(str(self.building_time_file_path),
                                str(self.complete_time_file_path))
                    if self.save_data is True:
                        self.post_blast_analysis(self.project)
                    # TODO-ROB Archive function here