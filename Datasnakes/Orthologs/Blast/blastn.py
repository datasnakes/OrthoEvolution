"""Optimized for use with local/standalone NCBI BLAST 2.6.0."""
import os
import shutil
import subprocess
import contextlib
from datetime import datetime as d
from pathlib import Path
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline

from Datasnakes.Manager import config
from Datasnakes.Orthologs.Blast.comparative_genetics_files import CompGenFiles
from Datasnakes.Orthologs.Blast.utils import gene_list_config, map_func
from Datasnakes.Tools.utils import makedirectory


# TODO-ROB: Find packages for script timing and analysis


class CompGenBLASTn(CompGenFiles):
    """Use CompGenBLASTn to search nucleotide databases using a nucleotide query.
    This class currently only works with the standalone blast.
    """

    def __init__(self, project, template=None, save_data=True, **kwargs):
        """Inherit from the BLASTing Template."""
        super().__init__(project=project, template=template, save_data=save_data, **kwargs)
        # # TODO-ROB Add taxon parameter
        # Manage Directories
        self.home = Path(os.getcwd())
        self.__gi_list_path = self.project_database / Path('gi_lists')

        makedirectory(self.__gi_list_path)

        # # Initialize Logging
        # self.__blastn_log = LogIt.blastn()
        #df = LogIt()
        self.date_format = '%a %b %d at %I:%M:%S %p %Y'
        # self.get_time = time.time  # To get the time use 'get_time()'
        # TODO-ROB:  Add a query organism variable
        self.query_gi_dict = {}
        self.removed_genes = []
        self.current_gene_list = []
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

    def blast_config(self, query_accessions, query_organism, auto_start=False):
        """Configure everything for a BLAST.
        First the accession file, and gene list is configured.
        """
        self.blastn_log.info(
            '***********************************BLAST CONFIG START************ \
            ***********************\n\n\n')
        self.blastn_log.info('Configuring the accession file...')

        # Update the gene_list based on the existence of a incomplete blast
        # file
        gene_list = gene_list_config(self.building_file_path, self.data, self.gene_list, self.taxon_dict, self.blastn_log)
        if gene_list is not None:
            start = len(self.blast_human) - len(gene_list)  # What gene to start blasting
            query_accessions = self.blast_human[start:]  # Reconfigured query
            # Reconfigure the gene_list to reflect the existing accession info
            self.current_gene_list = gene_list
        else:
            self.current_gene_list = self.gene_list

        # Create GI lists
        self.blastn_log.info("Configuring GI list using the taxonomy id and the blastdbcmd tool.")
        gi_list_config(self.__gi_list_path, self.research_path, self.taxon_ids, config)
        # Get GI (stdout) and query sequence (FASTA format)
        self.blastn_log.info("Generating directories.")
        self.blastn_log.info("Extracting query gi number to stdout and "
                             "query refseq sequence to a temp.fasta file from BLAST database.")
        # Iterate the query accessions numbers
        for query in query_accessions:
            gene = self.acc_dict[query][0][0]
            gene_path = self.raw_data / Path(gene) / Path('BLAST')
            # Create the directories for each gene
            try:
                Path.mkdir(gene_path, exist_ok=True, parents=True)
                self.blastn_log.info("Directory Created: %s" % gene)
                self.blastn_log.info("\n")
            except FileExistsError:
                self.blastn_log.info("Directory already exists: %s" % gene)

            # Save sequence data in FASTA file format and print the gi number to stdout with a custom BLAST extraction
            # https://www.ncbi.nlm.nih.gov/books/NBK279689/#_cookbook_Custom_data_extraction_and_form_
            # TODO-SDH Combine these BLAST extractions???
            fmt = {'query': query, 'temp fasta': str(gene_path / Path('temp.fasta'))}
            # Temporary fasta file created by the blastdbcmd
            fasta_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %f -out {temp fasta}".format(**fmt)
            fasta_status = subprocess.call([fasta_setup], shell=True)
            # Check the blast databases to see if the query accession even exists
            gi_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %g".format(**fmt)
            gi_status = subprocess.call([gi_setup], shell=True)
            # TODO-ROB:  Add function to add the gi numbers to the dataframe/csv-file,
            # TODO-ROB: and add a check function to see if thats already there
            # Check the status of the custom blast data extraction
            if gi_status == 0 or fasta_status == 0:  # One of the commands was successful.
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
                elif fasta_status != 0:
                    self.blastn_log.error(
                        "FASTA sequence for %s not found in the BLAST extraction" %
                        query)
                    self.blastn_log.error(
                        "Removing %s from the BLAST list..." % gene)
                    self.gene_list.remove(gene)
                    self.removed_genes.append(gene)
                    continue
                else:
                    pass  # Both commands were successful
            else:  # Both commands failed
                self.blastn_log.error(
                    "FASTA sequence and GI number for %s not found in the custom BLAST extraction." %
                    query)
                self.blastn_log.error(
                    "Removing %s from the BLAST list..." %
                    gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)
                with contextlib.suppress(ValueError):
                    self.current_gene_list.remove(gene)
                continue

        self.blastn_log.info('Configured gene list: %s\n\n\n' % self.current_gene_list)
        self.blastn_log.info(
            '************************************BLAST CONFIG END************************************\n\n\n')
        if auto_start is True:
            # Automatically begin BLASTING after the configuration
            self.blasting(
                genes=self.current_gene_list,
                query_organism=query_organism,
                pre_configured=auto_start)

    def blast_xml_parse(self, xml_path, gene, organism):
        """Parse the XML file created by the BLAST."""
        accession = gi = raw_bitscore = description = None
        record_dict = {}
        self.blastn_log.info("Parsing %s to find the best accession number." % xml_path)
        maximum = 0
        file_path = str(Path(xml_path))
        with open(file_path, 'r') as blast_xml:
            blast_qresult = SearchIO.read(blast_xml, 'blast-xml')
            mapped_qresult = blast_qresult.hit_map(map_func)  # Map the hits
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
            record_dict['gene'] = gene
            record_dict['organism'] = organism
            record_dict['accession'] = accession
            record_dict['gi'] = gi
            record_dict['raw bitscore'] = raw_bitscore
            record_dict['description'] = description
            self.blastn_log.info("Accession:  %s" % accession)
            self.blastn_log.info("GI number: {}".format(str(gi)))
            self.blastn_log.info("Raw bitscore: %s" % raw_bitscore)
            self.blastn_log.info("Description: %s" % description)
            self.add_accession(gene, organism, accession)

    def blasting(self, genes=None, query_organism=None, pre_configured=False):
        """Configure the BLAST."""
        if pre_configured is False:
            query = self.df[query_organism].tolist()
            self.blast_config(query_accessions=query, query_organism=query_organism, auto_start=True)
            genes = self.current_gene_list  # Gene list populated by blast_config.
        elif pre_configured is True:
            genes = genes

        self.blastn_log.info(
            "------------------------------------------------------------------")
        self.blastn_log.info("The script name is str(os.path.basename(__file__)).")
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
                gene_path = self.raw_data / Path(gene) / Path('BLAST')
                files = os.listdir(str(gene_path))
                xml = '%s_%s.xml' % (gene, organism)
                xml_path = gene_path / Path(xml)

                # Initialize configuration variables
                taxon_id = self.taxon_dict[organism]
                taxon_gi_file = str(taxon_id) + "gi"
                taxon_gi_path = self.__gi_list_path / Path(taxon_gi_file)

                if xml in files:
                    self.blast_xml_parse(xml_path, gene, organism)
                else:
                    self.blastn_log.warning("\n\n\n*******************BLAST START*******************")
                    start_time = self.get_time()
                    self.blastn_log.info("The start time is %s" % start_time)
                    self.blastn_log.info("The current gene is %s (%s)." % (gene, self.tier_dict[gene]))
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