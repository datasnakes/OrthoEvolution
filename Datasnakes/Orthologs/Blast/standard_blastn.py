"""Optimized for use with local/standalone NCBI BLAST 2.6.0 and higher"""
import os
import shutil
import subprocess
import contextlib
from datetime import datetime as d
from pathlib import Path
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline

from Datasnakes.Orthologs.Blast.comparative_genetics_files import CompGenFiles
from Datasnakes.Orthologs.Blast.utils import gene_list_config, map_func


class StandardBlastN(CompGenFiles):
    """Combines Project Management features with Blasting."""
    def __init__(self, template_file, **kwargs):
        """This class inherits from the CompGenFiles class.

        This class utilizes it's parent classes to search a standalone
        Blast database for specific orthologs of a gene using a query organism
        (usually human).  The best hits from the Blast are filtered for the
        best option in order to get the most accurate accession numbers for
        downstream analysis.

        :param kwargs:
        """
        super().__init__(project=None, template=template_file, save_data=True, **kwargs)

        # Create a date format
        self._fmt = '%a %b %d at %I:%M:%S %p %Y'
        self.date_format = str(d.now().strftime(self._fmt))

        # Ensure paths are set
        self.environment_vars = dict(os.environ)
        if 'BLASTDB' not in self.environment_vars.keys():
            msg = "BLASTDB is not set in your path."
            raise EnvironmentError(msg)
        elif 'WINDOW_MASKER_PATH' not in self.environment_vars.keys():
            msg = "WINDOW_MASKER_PATH is not set in your path."
            raise EnvironmentError(msg)

        # Manage Directories
        self.home = Path(os.getcwd())

        self.query_gi_dict = {}
        self.removed_genes = []
        self.current_gene_list = []

        self.blastn_log.info("These are the organisms: %s" % str(self.org_list))
        self.blastn_log.info("These are the genes: %s" % str(self.gene_list))
        self.blastn_log.info("These are the taxonomy ids: %s" % str(self.taxon_ids))

        # For completed blast files
        self.complete_file = self.project + '_MAF.csv'
        self.complete_file_path = self.data / Path(self.complete_file)
        self.complete_time_file = self.project + '_TIME.csv'
        self.complete_time_file_path = self.data / Path(self.complete_time_file)

    def blast_config(self, query_accessions, query_organism, auto_start=False):
        """This method configures everything for our BLAST workflow.

        It configures the accession file, which works with
        interrupted Blasts.  It configures a gene_list for blasting the right
        genes.

        :param query_accessions:  A list of query accession numbers.
                                  Each gene needs one from the same organism.
        :param query_organism:  The name of the query organism for post
                                configuration.
        :param auto_start:  A flag that determines whether the blast
                            starts automatically.
        :return:
        """
        sep = 20 * '*'
        self.blastn_log.info(sep + 'BLAST CONFIG START' + sep + '\n\n\n')
        self.blastn_log.info('Configuring the accession file.')

        # CONFIGURE and UPDATE the gene_list based on the existence of an
        # incomplete blast file
        gene_list = gene_list_config(self.building_file_path, self.data,
                                     self.gene_list, self.taxon_dict,
                                     self.blastn_log)
        if gene_list is not None:
            # What gene to start blasting
            start = len(self.blast_human) - len(gene_list)

            # Reconfigured query
            query_accessions = self.blast_human[start:]

            # Reconfigure the gene_list to reflect the existing accession info
            self.current_gene_list = gene_list
        else:
            self.current_gene_list = self.gene_list

        # Get GI (stdout) and query sequence (FASTA format)
        self.blastn_log.info("Creating directories.")
        self.blastn_log.info("Extracting query refseq sequence to a temp.fasta file from BLAST database.")

        # Iterate the query accessions numbers
        for query in query_accessions:
            gene = self.acc_dict[query][0][0]
            gene_path = self.raw_data / Path(gene) / Path('BLAST')
            # Create the directories for each gene
            try:
                Path.mkdir(gene_path, exist_ok=True, parents=True)
                self.blastn_log.info("Directory created: %s" % gene)
            except FileExistsError:
                self.blastn_log.info("Directory exists: %s" % gene)

            # Save sequence data in FASTA file format and print the gi number
            # to stdout with a custom BLAST extraction
            # https://www.ncbi.nlm.nih.gov/books/NBK279689/#_cookbook_Custom_data_extraction_and_form_
            # TODO-SDH Combine these BLAST extractions???
            fmt = {'query': query, 'temp fasta': str(gene_path / Path('temp.fasta'))}
            # Temporary fasta file created by the blastdbcmd
            fasta_setup = "blastdbcmd -entry {query} -db refseq_rna -outfmt %f -out {temp fasta}".format(**fmt)
            fasta_status = subprocess.call([fasta_setup], shell=True)

            # Check the status of the custom blast data extraction
            if fasta_status != 0:  # Command was successful.
                self.blastn_log.error("FASTA sequence for %s not found in the BLAST extraction" % query)
                self.blastn_log.error("Removing %s from the BLAST list." % gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)

            else:  # Command failed
                self.blastn_log.error("FASTA sequence not found in the custom BLAST extraction." % query)
                self.blastn_log.error("Removing %s from the BLAST list." % gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)
                with contextlib.suppress(ValueError):
                    self.current_gene_list.remove(gene)
                continue

        self.blastn_log.info('Configured gene list: %s' % self.current_gene_list)
        self.blastn_log.info(sep + 'BLAST CONFIG END' + sep)
        if auto_start is True:
            # Automatically begin BLASTING after the configuration
            self.blasting(genes=self.current_gene_list,
                          query_organism=query_organism,
                          pre_configured=auto_start)

    def blast_xml_parse(self, xml_path, gene, organism):
        """Parse the blast XML record get the best hit accession number.

        :param xml_path:  Absolute path to the blast record.
        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :return:  Returns one accession number in the building accession file.
        """
        accession = gi = raw_bitscore = description = None
        record_dict = {}
        xmlsplit = xml_path.split('/')
        curxmlfile = xmlsplit[-1]
        self.blastn_log.info("Parsing %s to find the best accession number." % curxmlfile)
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
                        # TODO Add another check here???
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
            self.blastn_log.info("Finished parsing the Blast XML file.")
            self.blastn_log.info("The best accession has been selected.")
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
        """Run NCBI's blastn.

        This method actually performs NCBI's blastn.
        It requires configuring before it can be utilized.

        :param genes:  Gene of interest.
        :param query_organism:  Query organism.
        :param pre_configured:  Determines if the blast needs configuring.
        :return:
        """
        ast = 10 * '*'  # Asterisk separator
        if pre_configured is False:
            query = self.df[query_organism].tolist()
            self.blast_config(query_accessions=query,
                              query_organism=query_organism,
                              auto_start=True)
            # Gene list populated by blast_config
            genes = self.current_gene_list
        elif pre_configured is True:
            genes = genes

        self.blastn_log.info('The script began on {}'.format(str(d.now().strftime(self.date_format))))

        # TIP Steps to a bulk blast
        # 1.  Iterate the gene of interest
        # 2.  Iterate the organisms of interest
        # 3.  Perform BlastN
        for gene in genes:
            for organism in self.org_list:
                # Skip the query organism
                if organism == query_organism:
                    continue
                # Initialize output variables
                gene_path = self.raw_data / Path(gene) / Path('BLAST')
                files = os.listdir(str(gene_path))
                xml = '%s_%s.xml' % (gene, organism)
                xml_path = str(gene_path / Path(xml))

                # Initialize configuration variables
                taxon_id = self.taxon_dict[organism]


                if xml in files:
                    self.blast_xml_parse(xml_path, gene, organism)
                else:
                    self.blastn_log.warning(ast + "BLAST START" + ast)
                    start_time = self.get_time()
                    self.blastn_log.info("The start time is %s" % start_time)
                    self.blastn_log.info("The current gene is %s (%s)." % (gene, self.tier_dict[gene]))
                    self.blastn_log.info("The current organisms is %s (%s)." % (organism, taxon_id))

                    with open(xml_path, 'w') as blast_xml:
                        # Set up blast parameters
                        query_seq_path = str(gene_path / Path('temp.fasta'))

                        # Use Biopython's NCBIBlastnCommandline tool
                        result_handle = NcbiblastnCommandline(query=query_seq_path,
                                                              db="refseq_rna",
                                                              strand="plus",
                                                              evalue=0.001,
                                                              outfmt=5,
                                                              window_masker_taxid=taxon_id,
                                                              max_target_seqs=10,
                                                              task="blastn")

                        # Capture blast data
                        stdout_str, stderr_str = result_handle()
                        self.blastn_log.error(stderr_str)
                        blast_xml.write(stdout_str)
                        end_time = self.get_time()
                        elapsed_time = end_time - start_time
                        self.blastn_log.info("%s was create." % blast_xml.name)
                        self.blastn_log.info("The end time is %s." % end_time)
                        self.blastn_log.info("The BLAST took %s." % elapsed_time)
                    self.blastn_log.warning(ast + "BLAST END" + ast)
                    self.add_blast_time(gene, organism, start_time, end_time)
                    self.blast_xml_parse(xml_path, gene, organism)

                # If the BLAST has gone through all orthologs then create a
                # master accession file.
                if gene == genes[-1] and organism == self.org_list[-1]:
                    shutil.copy(str(self.building_file_path),
                                str(self.complete_file_path))
                    shutil.copy(str(self.building_time_file_path),
                                str(self.complete_time_file_path))
                    if self.save_data is True:
                        self.post_blast_analysis(self.project)
                        self.blastn_log.info("Post blast analysis is complete")

        self.blastn_log.info("BLAST has completed. Check your output data.")
