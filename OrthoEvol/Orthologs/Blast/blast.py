"""Optimized for use with local/standalone NCBI BLAST 2.8.1"""
import os
import shutil
import contextlib
from subprocess import run, PIPE, CalledProcessError
from datetime import datetime as d
import time
from pathlib import Path
import logging

from Bio.Application import ApplicationError
from Bio import SearchIO  # Used for parsing and sorting XML files.

from OrthoEvol.Orthologs.Blast.blastn_wrapper import NcbiblastnCommandline
from OrthoEvol.Orthologs.Blast.comparative_genetics import ComparativeGenetics
# Other
from xml.etree.ElementTree import ParseError


class BaseBlastN(ComparativeGenetics):
    """Base BlastN class."""

    def __init__(self, project, method, acc_file, template=None, save_data=True,
                 verbose=False, **kwargs):
        """This class inherits from the ComparativeGenetics class.
        This class utilizes it's parent classes to search a standalone
        Blast database for specific orthologs of a gene using a query organism
        (usually human).  The best hits from the Blast are filtered for the
        best option in order to get the most accuarate accession numbers for
        downstream analysis.

        :param project:  The project name.
        :param method:  Method used for blasting. (1, 2, or None)
        :param template:  The accession file template.
        :param save_data:  A flag for saving the post_blast data to an excel file.
        :param quiet:  A flag for determining the level of logging verbosity.
        :param kwargs:"""

        super().__init__(project=project, method=method, acc_file=acc_file,
                         template=template, save_data=save_data, verbose=verbose, **kwargs)

        if verbose:
            self.blastn_log.setLevel(logging.DEBUG)
        else:
            self.blastn_log.setLevel(logging.INFO)

        # Method variable
        self.method = method
        self.acc_file = acc_file

        # Create a date format
        self._fmt = '%a %b %d at %I:%M:%S %p %Y'
        self.date_format = str(d.now().strftime(self._fmt))

        # Ensure paths are set
        self.environment_vars = dict(os.environ)
        if 'BLASTDB' not in self.environment_vars.keys():
            msg = "BLASTDB is not set in your path."
            raise EnvironmentError(msg)

        # Manage Directories
        self.home = Path(os.getcwd())
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

        self.blastn_parameters, self.query_config = BaseBlastN.select_method(method=self.method)

    def _make_blast_dir(self, gene, path):
        """Create a blast directory for a gene.

        :param gene:  A gene that will be run with blastn.
        :param path:  The path of the gene directory named with the gene name.
        """
        try:
            Path.mkdir(path, exist_ok=True, parents=True)
            self.blastn_log.debug("Directory created for %s" % gene)
        except FileExistsError:
            self.blastn_log.debug("Directory exists for %s" % gene)

    def _create_temp_fasta(self, query, gene, query_config):
        """Create a temporary fasta file using blastdbcmd.

        :param query: The blast query or reference accession.
        :param gene: The gene name of the accession.
        :param query_config: A configuration dict for the blastdbcmd string.
        """
        try:
            blastdbcmd_query = "blastdbcmd -entry {query} -db {db} -outfmt %f -out {temp fasta}".format(
                **query_config)
            blastdbcmd_status = run(blastdbcmd_query, stdout=PIPE,
                                    stderr=PIPE, shell=True, check=True)
            self.blastn_log.info(
                "Extracted query refseq sequence to a temp.fasta file from BLAST database.")
        except CalledProcessError as err:
            self.blastn_log.error(err.stderr.decode('utf-8'))
        else:
            if blastdbcmd_status.returncode != 0:  # Command was not successful.
                self.blastn_log.error(
                    "FASTA sequence for %s not found in the BLAST extraction" % query)
                self.blastn_log.error("Removing %s from the BLAST list." % gene)
                self.gene_list.remove(gene)
                self.removed_genes.append(gene)

                with contextlib.suppress(ValueError):
                    self.current_gene_list.remove(gene)

    @staticmethod
    def select_method(method=1):
        """Select a method for running blastn.

        :param method: The blast method to use. Either 1, 2, or None.
        """
        # Local blast using taxonomy ids
        # TIP: This is the fastest version and recommended for blasting with
        # multiple species and multiple genes.
        if method == 1:
            blastn_parameters = {'query': '', 'db': 'refseq_rna_v5',
                                 'taxids': '', 'strand': 'plus',
                                 'evalue': 0.01, 'outfmt': 5,
                                 'max_target_seqs': 10, 'task': 'blastn'}
            query_config = {'query': '',
                            'db': 'refseq_rna_v5',
                            'temp fasta': ''}
        # Remote blast with entrez_query
        # XXX: This cannot be used with taxids.
        elif method == 2:
            blastn_parameters = {'query': '', 'entrez_query': '',
                                 'db': 'refseq_rna',
                                 'strand': 'plus', 'evalue': 0.01,
                                 'outfmt': 5, 'max_target_seqs': 10,
                                 'task': 'blastn', 'remote': 'True'}
            query_config = {'query': '',
                            'db': 'refseq_rna',
                            'temp fasta': ''}
        # Default local blast
        # XXX Unless you are have a very simple blast....
        # i.e. only a query fasta
        elif method is None or method == "":
            blastn_parameters = {'query': '', 'db': 'refseq_rna',
                                 'strand': 'plus', 'evalue': 0.01,
                                 'outfmt': 5, 'max_target_seqs': 500,
                                 'task': 'blastn'}
            query_config = {'query': '',
                            'db': 'refseq_rna',
                            'temp fasta': ''}
        else:
            raise ValueError('%s is not a blast method.' % method)
        return blastn_parameters, query_config

    def blastn_wrapper(self, gene, organism, parameters, xml_path, gene_path):
        """Run the Ncbiblastncommandline wrapper modified by this package."""
        try:
            self.blastn_log.info('Blast run has started.')
            start_time = self.get_time()
            with open(xml_path, 'w') as blast_xml:
                # Update blastn parameters
                self.blastn_parameters.update(parameters)

                # Use Biopython's NCBIBlastnCommandline tool
                result_handle = NcbiblastnCommandline(**self.blastn_parameters)
                # Capture the standard output
                stdout_str, _ = result_handle()

                # Write the stdout_str to the xml file
                blast_xml.write(stdout_str)
                end_time = self.get_time()
                elapsed_time = round(end_time - start_time, 2)
                self.blastn_log.debug("%s was created." % blast_xml.name)
                self.blastn_log.info("The BLAST took %ss." % elapsed_time)
                self.blastn_log.debug('Blast run has ended.')
                self.add_blast_time(gene, organism, start_time, end_time)
                blast_xml.close()

                # Parse the blast output file
                self.parse_xml(xml_path, gene, organism)

                # Make blast sleep if server blast
                if self.method == 2:
                    time.sleep(30)  # Prevents server lockout

        # Catch either ApplicationError or ParseError
        except (ApplicationError, ParseError) as err:
            try:
                errmsg = 'Parse Error - %s' % err.msg
            except AttributeError:
                errmsg = err.stderr
            self.blastn_log.error(errmsg)
            os.remove(xml_path)
            self.blastn_log.info('%s was deleted.' % xml_path)

        # Catch either KeyboardInterrupt or SystemExit
        except (KeyboardInterrupt, SystemExit):
            msg = "Keyboard interrupt or System Exit."
            self.blastn_log.error(msg)
            os.remove(xml_path)
            self.blastn_log.error('%s was deleted.' % xml_path)
            raise

    def configure(self, query_accessions, query_organism, auto_start=False):
        """This method configures everything for our BLAST workflow.
        It configures the accession file, which works with interrupted Blasts.
        It configures a gene_list for blasting the right genes.

        :param query_accessions:  A list of query accession numbers.  Each gene
                                  needs one from the same organism.
        :param query_organism:  The name of the query organism for post
                                configuration.
        :param auto_start:  A flag that determines whether the blast starts
                            automatically. (Default value = False)
        :return:
        """

        self.blastn_log.debug('Blast configuration has begun.')
        self.blastn_log.debug('Configuring the accession file.')

        # CONFIGURE and UPDATE the gene_list based on the existence of an
        # incomplete blast file
        gene_list = self.blast_utils.gene_list_config(self.building_file_path, self.data,
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

        # Iterate through the query accessions numbers
        for query in query_accessions:
            gene = self.acc_dict[query][0][0]
            gene_path = self.raw_data / Path(gene) / Path('BLAST')

            # Create the directories for each gene
            self._make_blast_dir(gene=gene, path=gene_path)

            # Save sequence data in FASTA file format and print the gi number
            # to stdout with a custom BLAST extraction
            # https://www.ncbi.nlm.nih.gov/books/NBK279689/#_cookbook_Custom_data_extraction_and_form_
            self.query_config.update({'query': query,
                                      'temp fasta': str(gene_path / Path('temp.fasta'))})

            # Create a temporary fasta file using blastdbcmd
            self._create_temp_fasta(query=query, gene=gene,
                                    query_config=self.query_config)

        self.blastn_log.debug('Configured gene list: %s' % self.current_gene_list)
        self.blastn_log.debug('Blast configuration has ended.')
        if auto_start:  # Automatically run blast after the configuration
            self.runblast(genes=self.current_gene_list,
                          query_organism=query_organism,
                          pre_configured=auto_start)

    def create_maf(self):
        """Create a master accessions file if the blast has completed."""
        try:
            shutil.copyfile(str(self.building_file_path),
                            str(self.complete_file_path))
            shutil.copyfile(str(self.building_time_file_path),
                            str(self.complete_time_file_path))
        except FileNotFoundError as fnfe:
            msg = "Your blast did not create building time files."
            self.blastn_log.error(fnfe)
            self.blastn_log.error(msg)
            raise BlastFailure("Blast has failed to generate "
                               "data. Check logs for errors.")

    def parse_xml(self, xml_path, gene, organism):
        """Parse the blast XML record get the best hit accession number.

        :param xml_path:  Absolute path to the blast record.
        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :return:  Returns one accession number in the building accession file.
        """
        accession = gi = raw_bitscore = description = None
        record_dict = {}
        xmlsplit = xml_path.split(os.sep)
        curxmlfile = xmlsplit[-1]
        self.blastn_log.info("Parsing %s to find the best accession number." % curxmlfile)
        maximum = 0
        file_path = str(Path(xml_path))

        # Read the blast_xml output file
        with open(file_path, 'r') as blast_xml:
            blast_qresult = SearchIO.read(blast_xml, 'blast-xml')
            # Map the hits
            mapped_qresult = blast_qresult.hit_map(self.blast_utils.map_func)  # Map the hits
            # Extract hit information from xml file.
            for hit in mapped_qresult:
                for hsp in hit.hsps:
                    # Find the highest scoring hit for each gene
                    if hsp.bitscore_raw > maximum:
                        # If the gene is a predicted non-coding RefSeq gene
                        # then go the the next hit
                        # https://en.wikipedia.org/wiki/RefSeq
                        if "xr" in str(hit.id.lower()):
                            self.blastn_log.info("Encountered a predicted(X*_00000) "
                                                 "non-coding (XR_000000)(ncRNA) RefSeq "
                                                 "gene.  Moving to the next hit.")
                            break
                        else:
                            # If the gene is acceptable then add it to the gene list
                            # Lowercase means LOC1223214 is the name of the gene
                            # TODO:  Change this check
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
            blast_xml.close()
            self.blastn_log.info("Finished parsing the Blast XML file.")
            self.blastn_log.info("The best accession has been selected.")

            # Add xml parsed xml data to record dictionary
            record_dict['gene'] = gene
            record_dict['organism'] = organism
            record_dict['accession'] = accession
            record_dict['gi'] = gi
            record_dict['raw bitscore'] = raw_bitscore
            record_dict['description'] = description

            # Log the parsed xml data
            self.blastn_log.info("Accession:  %s" % accession)
            self.blastn_log.info("GI number: {}".format(str(gi)))
            self.blastn_log.info("Raw bitscore: %s" % raw_bitscore)
            self.blastn_log.info("Description: %s" % description)
            self.add_accession(gene, organism, accession)

    def runblast(self, genes=None, query_organism=None, pre_configured=False):
        """Run NCBI's blastn.
        This method actually performs NCBI's blastn.
        It requires configuring before it can be utilized.

        :param genes:  Gene of interest. (Default value = None)
        :param query_organism:  Query organism. (Default value = None)
        :param pre_configured:  Determines if the blast needs configuring. (Default value = False)
        :return:
        """
        if not pre_configured:
            query = self.df[query_organism].tolist()
            self.configure(query_accessions=query,
                           query_organism=query_organism,
                           auto_start=True)
            # Gene list created by blast_config.
            genes = self.current_gene_list

        if len(genes) < 1:
            self.blastn_log.fatal('There are no genes to run with blastn.')

        else:
            self.blastn_log.debug('The blast began on {}'.format(self.date_format))

            # TIP Steps to run blastn for orthology inference
            # 1.  Iterate over genes of interest
            # 2.  Iterate over organisms of interest
            # 3.  Perform blastn
            # 4.  Parse blastn output for best hit accession
            for gene in genes:
                for organism in self.org_list:
                    # Skip the query organism
                    if organism == query_organism:
                        continue
                    # Initialize output variables
                    gene_path = self.raw_data / Path(gene) / Path('BLAST')
                    xml = '%s_%s.xml' % (gene, organism)
                    xml_path = str(gene_path / Path(xml))

                    # Initialize configuration variables
                    taxon_id = self.taxon_dict[organism]

                    if xml in os.listdir(str(gene_path)):
                        # Try to parse your xml file. If it is not parsable,
                        # catch the error.
                        try:
                            self.parse_xml(xml_path, gene, organism)
                        # The error will likely be the result of an empty file.
                        except ParseError as err:
                            parse_error_msg = 'Parse Error - %s' % err.msg
                            self.blastn_log.error(parse_error_msg)
                            os.remove(xml_path)
                            self.blastn_log.error('%s was deleted' % xml)

                    else:
                         # Set up blast parameters
                        query_seq_path = str(gene_path / Path('temp.fasta'))
                        # Add blastn parameters for each method to dict
                        if self.method == 1:
                            update_dict = {'query': query_seq_path,
                                           'taxids': taxon_id}
                        elif self.method == 2:
                            update_dict = {'query': query_seq_path,
                                           'entrez_query': organism}
                        else:
                            update_dict = {'query': query_seq_path}

                        self.blastn_log.debug('Blast run has started.')
                        self.blastn_log.info("Current gene: %s (%s)." %
                                             (gene, self.tier_dict[gene]))
                        self.blastn_log.info("Current organism: %s (%s)." % (organism, taxon_id))
                        # Run Ncbiblastncommandline
                        self.blastn_wrapper(xml_path=xml_path, gene_path=gene_path,
                                            parameters=update_dict, gene=gene,
                                            organism=organism)

                    # If the BLAST has gone through all orthologs then create a
                    # master accession file.
                    if gene == genes[-1] and organism == self.org_list[-1]:
                        self.create_maf()
                        if self.save_data:
                            self.post_blast_analysis(self.project)
                            self.blastn_log.info("Post-blast analysis is complete")

                        # TODO Archive function here
            self.blastn_log.info("Blast has completed. Check your output.")


class OrthoBlastN(BaseBlastN):
    """OrthoBlastN provides a preconfigured and orthology optimized version
    of the BaseBlastN class."""

    def __init__(self, project="orthology-gpcr", method=1, template=None,
                 save_data=True, acc_file="gpcr.csv", copy_from_package=True,
                 **kwargs):
        """This class inherits from the BaseBlastN class.

        This class utilizes it's parent classes to search a standalone
        Blast database for specific orthologs of a gene using a query organism
        (usually human).  The best hits from the Blast are filtered for the
        best option in order to get the most accuarate accession numbers for
        downstream analysis.

        :param project:  The project name (Default: 'orthology-gpcr')
        :param method: Method used for blasting. (Default: 1)
        :param template:  The accession file template.
        :param save_data:  A flag for saving the post_blast data to an excel file.
        :param acc_file: The accession file to use. (Default: 'gpcr.csv')
        :param copy_from_package: Copy the acc_file from the package. (Default: True)
        :param kwargs:
        """
        # Set values for methods to prevent using a config.
        self.taxon_file = None
        self.__post_blast = True
        self.go_list = None
        self.project_path = None
        self.proj_mana = None
        self.acc_file = self.MAF = acc_file
        self.copy_from_package = copy_from_package

        # Initialize class
        super().__init__(project=project, method=method, template=template,
                         save_data=save_data, acc_file=self.acc_file,
                         copy_from_package=self.copy_from_package,
                         MAF=self.MAF,
                         taxon_file=self.taxon_file,
                         post_blast=self.__post_blast,
                         go_list=self.go_list, project_path=self.project_path,
                         proj_mana=self.proj_mana, **kwargs)

    def run(self):
        """Run the blast using a default configuration."""
        self.configure(self.blast_human, self.species, auto_start=True)


class BlastFailure(BaseException):
    """Simple Blast Exception Class."""
    pass
