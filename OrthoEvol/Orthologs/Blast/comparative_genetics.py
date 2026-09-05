"""Comparative Genetics"""
# Standard Library
import copy
import os
import random
import shutil
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

# Other
import pandas as pd
from ete3 import NCBITaxa

# OrthoEvol
from OrthoEvol.Manager.config import data
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.resources import package_resource_path
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.utilities import FullUtilities

# from pandas import ExcelWriter
# NCBITaxa().update_taxonomy_database()

# TODO: Create function for archiving and multiple runs (this can go
# into the Management class)

HGNC_COMPLETE_SET_URL = (
    "https://storage.googleapis.com/public-download-files/"
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
)
HGNC_FIELDS = (
    "hgnc_id",
    "symbol",
    "name",
    "status",
    "locus_type",
    "entrez_id",
    "ensembl_gene_id",
    "refseq_accession",
)


def _count_worksheet(
    counts: Mapping[str, int | float],
) -> pd.DataFrame | None:
    """Create a consistently labeled count table when counts are available."""
    if not counts:
        return None
    return pd.DataFrame.from_dict(counts, orient="index", columns=["Count"])


def _mapping_worksheet(
    values: Mapping[str, object],
) -> pd.DataFrame | None:
    """Create a worksheet only when the corresponding analysis has results."""
    if not values:
        return None
    return pd.DataFrame.from_dict(values, orient="index")


def _duplicate_group_worksheet(
    groups: Mapping[str, Mapping[str, Sequence[str]]],
) -> pd.DataFrame | None:
    """Collect duplicate groups while retaining their outer entity labels."""
    if not groups:
        return None
    grouped_values = {
        entity: list(accession_groups.values())
        for entity, accession_groups in groups.items()
    }
    return pd.DataFrame.from_dict(grouped_values, orient="index").T


def _missing_worksheets(
    missing_records: Mapping[str, Mapping[str, object]],
    details_key: str,
    count_sheet_name: str,
    details_sheet_name: str,
) -> dict[str, pd.DataFrame]:
    """Separate missing-item details from their per-entity counts."""
    if not missing_records:
        return {}

    required_keys = {"count", details_key}
    for entity, record in missing_records.items():
        missing_keys = required_keys.difference(record)
        if missing_keys:
            missing_key_list = ", ".join(sorted(missing_keys))
            raise ValueError(
                f"Missing-data record {entity!r} requires: {missing_key_list}."
            )

    counts = {
        entity: record["count"] for entity, record in missing_records.items()
    }
    details = {
        entity: record[details_key]
        for entity, record in missing_records.items()
    }
    return {
        count_sheet_name: pd.DataFrame.from_dict(
            counts,
            orient="index",
            columns=["Count"],
        ),
        details_sheet_name: pd.DataFrame.from_dict(details, orient="index"),
    }


class BaseComparativeGenetics(object):
    """Base class in the Blast module."""

    __acc_filename = ''
    __paml_filename = ''
    __acc_path = ''
    __data = ''

    # Initialize Logging
    blastn_log = LogIt().default(logname="blastn", logfile=None)

    # TODO:  CREATE PRE-BLAST and POST-BLAST functions
    def __init__(
        self,
        project: str | None = None,
        project_path: str | Path | None = os.getcwd(),
        acc_file: str | None = None,
        taxon_file: str | Path | None = None,
        ref_species: str | None = None,
        pre_blast: bool = False,
        post_blast: bool = True,
        hgnc: bool | str | Path = False,
        proj_mana: ProjectManagement | dict[str, Any] | None = None,
        copy_from_package: bool = False,
        **kwargs: Any,
    ) -> None:
        """This is the base class for the Blast module.

        It parses an accession file in order to provide easy handling for data.

        The .csv accession file contains the following header info:
            * "Tier" - User defined.
            * "Gene" - HUGO Gene Nomenclature Committee(HGNC) symbol for the genes of interest.
            * Query Organism - A well annotated query organism.
            * Other organisms - The other headers are Genus_species of other taxa.

        The organisms are taken from:
        ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/
        The genes are taken from:
        http://www.guidetopharmacology.org/targets.jsp.
        The API gives the user access to their data in a higher level for
        downstream processing or for basic observation of the data.

        :param project:  The name of the project.
        :param project_path:  The location of the project, which is generally
                            defined by the ProjectManagement configuration.
        :param acc_file:  The name of the accession file.
        :param taxon_file:  A file that contains an ordered list of taxonomy ids.
        :param ref_species: A reference species or organism for the blast query.
        :param pre_blast:  A flag that gives the user access to an API that
                        contains extra information about their genes using the
                        mygene package.
        :param post_blast:  A flag that is used to handle a BLAST result file,
                            which returns information about misssing
                            data, duplicates, etc.
        :param hgnc:  Enable HGNC annotation with the current complete dataset,
                      or provide a local HGNC TSV path or alternate URL.
        :param proj_mana:  This parameter is used to compose (vs inherit) the
                           ProjectManagement class with the ComparativeGenetics class.
                           This parameter allows the various blast classes to function with or
                           without the Manager module.
        :param copy_from_package: Copy a packaged accession file into the project.
        :param kwargs:  The kwargs here are generally used for standalone blasting or for development.
        :returns:  A pandas data-frame, pivot-table, and associated lists and dictionaries.
        """

        self.__pre_blast = pre_blast
        self.__post_blast = post_blast
        self.__hgnc_source = self._resolve_hgnc_source(hgnc)
        self.acc_file = acc_file
        self.project_path = project_path
        self.project = project
        self.ref_species = ref_species
        self.taxon_file = taxon_file
        self.proj_mana = proj_mana
        self.get_time = time.time
        self.sep = 50 * "*"
        self.blast_utils = FullUtilities()

        project, project_path = self._resolve_project_context(
            project,
            project_path,
        )

        self.blastn_log.debug('Project name: %s' % self.project)
        self.blastn_log.debug('Project path: %s' % self.project_path)

        self._apply_project_configuration(project, project_path, proj_mana)
        self._initialize_taxon_path()
        if self.acc_file is not None:
            self._copy_accession_file(copy_from_package)
        self.acc_filename = self.acc_file
        self._initialize_accession_data(acc_file)

    @staticmethod
    def _resolve_hgnc_source(
        hgnc: bool | str | Path,
    ) -> str | Path | None:
        """Resolve the opt-in flag without changing custom HGNC sources."""
        if hgnc is True:
            return HGNC_COMPLETE_SET_URL
        if hgnc:
            return hgnc
        return None

    def _resolve_project_context(
        self,
        project: str | None,
        project_path: str | Path | None,
    ) -> tuple[str | None, str | Path | None]:
        """Resolve instance paths and the inputs used by project composition."""
        if project_path and project:
            self.project_path = Path(project_path) / project
        elif project and not project_path:
            self.project_path = project
        elif not project and not project_path:
            # Retain the legacy generated name for unnamed standalone projects.
            four_ints = random.sample(range(1, 9), 4)
            self.project = "orthoevol" + "".join(str(value) for value in four_ints)
            self.project_path = os.getcwd()
        elif project_path and not project:
            existing_project_path = Path(project_path)
            self.project = existing_project_path.name
            project = self.project
            project_path = existing_project_path.parent
        return project, project_path

    def _apply_project_configuration(
        self,
        project: str | None,
        project_path: str | Path | None,
        proj_mana: ProjectManagement | dict[str, Any] | None,
    ) -> None:
        """Copy composed project attributes onto this analysis instance."""
        configured_self = self.blast_utils.attribute_config(
            cls=self,
            composer=proj_mana,
            checker=ProjectManagement,
            project=project,
            project_path=project_path,
        )
        for variable, attribute in configured_self.__dict__.items():
            setattr(self, variable, attribute)

    def _initialize_taxon_path(self) -> None:
        """Create the optional taxonomy path after project composition."""
        if self.taxon_file is not None:
            self.taxon_path = self.project_index / Path(self.taxon_file)

    def _copy_accession_file(self, copy_from_package: bool) -> None:
        """Copy the selected accession source into the project index."""
        if copy_from_package:
            accession_source = package_resource_path(data, self.acc_file)
        else:
            accession_source = self.acc_file
        shutil.copy(accession_source, str(self.project_index))

    def _initialize_accession_data(self, original_acc_file: str | None) -> None:
        """Initialize accession-backed state after the source file is copied."""
        if self.acc_file is None:
            self._initialize_empty_output_filenames()
            return

        self._initialize_accession_paths(original_acc_file)
        self._initialize_accession_collections()
        self._initialize_accession_frames()
        if self.__post_blast:
            self._initialize_post_blast_state()
        self._initialize_accession_views()

    def _initialize_accession_paths(
        self,
        original_acc_file: str | None,
    ) -> None:
        """Derive all accession-related input and output paths."""
        self.acc_sqlite_filename = Path(original_acc_file).stem + ".sqlite"
        self.acc_sqlite_tablename = Path(original_acc_file).stem.replace(".", "_")
        self.acc_csv_path = self.project_index / Path(self.acc_file)
        self.acc_sqlite_path = self.project_index / self.acc_sqlite_filename
        self.building_filename = f"{self.acc_file[:-4]}building.csv"
        self.building_file_path = self.data / self.building_filename
        self.building_time_filename = self.building_filename.replace(
            "building.csv",
            "building_time.csv",
        )
        self.building_time_file_path = self.data / self.building_time_filename
        self.mygene_filename = f"{self.project}_mygene.csv"
        self.mygene_path = self.data / self.mygene_filename
        self.hgnc_filename = f"{self.project}_hgnc.csv"
        self.hgnc_path = self.data / self.hgnc_filename

    def _initialize_accession_collections(self) -> None:
        """Create the mutable collections populated during BLAST analysis."""
        self.org_list = []
        self.ncbi_orgs = []
        self.org_count = 0
        self.taxon_ids = []
        self.taxon_orgs = []
        self.taxon_dict = {}
        self.gene_list = []
        self.gene_count = 0
        self.tier_list = []
        self.tier_dict = {}
        self.tier_frame_dict = {}
        self.acc_dict = {}
        self.acc_list = []
        self.blast_human = []
        self.blast_rhesus = []

    def _initialize_accession_frames(self) -> None:
        """Load the accession table and create its mutable working frames."""
        self.raw_acc_data = self.blast_utils.accession_sqlite2pandas(
            self.acc_sqlite_tablename,
            self.acc_sqlite_filename,
            path=self.project_index,
            acc_file=self.acc_file,
        )
        self.mygene_df = pd.DataFrame()
        self.hgnc_df = pd.DataFrame()
        self.header = self.raw_acc_data.axes[1].tolist()

        self.building = copy.deepcopy(self.raw_acc_data)
        del self.building["Tier"]
        del self.building[self.ref_species]
        self.building = self.building.set_index("Gene")

        self.building_time = copy.deepcopy(self.raw_acc_data)
        del self.building_time["Tier"]
        del self.building_time[self.ref_species]
        self.building_time = self.building_time.set_index("Gene")

    def _initialize_post_blast_state(self) -> None:
        """Create result containers only when post-BLAST analysis is enabled."""
        self.missing_dict = {}
        self.missing_genes = {}
        self.missing_organsims = {}
        self.missing_gene_count = 0
        self.missing_organsims_count = 0
        self.duplicated_dict = {}
        self.duplicated_accessions = {}
        self.dup_acc_count = {}
        self.duplicated_genes = {}
        self.dup_gene_count = {}
        self.duplicated_organisms = {}
        self.dup_org_count = {}
        self.duplicated_random = {}
        self.duplicated_other = {}
        self.time_dict = {}

    def _initialize_accession_views(self) -> None:
        """Create indexed views and populate the public lookup collections."""
        self.__data = self.raw_acc_data.set_index("Gene")
        self.df = self.__data
        self.pt = pd.pivot_table(
            copy.deepcopy(self.raw_acc_data),
            index=["Tier", "Gene"],
            aggfunc="first",
        )
        organism_columns = self.pt.axes[1].tolist()
        self.pt.columns = pd.Index(organism_columns, name="Organism")
        self.org_dict = self.df.loc[:, self.ref_species:].to_dict()
        self.gene_dict = self.df.T.to_dict()
        self.get_master_lists(self.__data)

    def _initialize_empty_output_filenames(self) -> None:
        """Retain output names when no accession dataset is configured."""
        self.building_filename = f"{self.project}_building.csv"
        self.building_time_filename = f"{self.project}_building_time.csv"

    @staticmethod
    def get_hgnc_gene_info(
        gene_symbols: Sequence[str],
        source: str | Path = HGNC_COMPLETE_SET_URL,
    ) -> pd.DataFrame:
        """Return HGNC records in the same order as the requested symbols."""
        requested_symbols = pd.Index(
            (symbol.strip().upper() for symbol in gene_symbols),
            dtype="string",
            name="query_symbol",
        )
        if requested_symbols.empty:
            return pd.DataFrame(columns=("query_symbol", *HGNC_FIELDS))

        hgnc_records = pd.read_csv(
            source,
            sep="\t",
            usecols=HGNC_FIELDS,
            dtype="string",
        ).set_index("symbol", drop=False)

        matched_records = hgnc_records.reindex(requested_symbols).reset_index(
            drop=True
        )
        matched_records.insert(0, "query_symbol", requested_symbols)
        return matched_records

    @staticmethod
    def get_file_list(file):
        """Turn csv column to list.

        :param file: Name of csv file.
        """
        file_data = pd.read_csv(file, header=None)
        file_list = list(file_data[0])
        return file_list

    def get_master_lists(self, df, csv_file=None):
        """Populate the organism and gene lists with a data frame.

        It will also populate pre-blast attributes (mygene) and post-blast
        attributes (missing and duplicates) under the proper conditions.

        :param df: The preferred way of utilizing the function is with a data-frame.
        :param csv_file: If a csv_file is given, then a data-frame will be
                         created by reinitializing the object.
                         (Default value = None)
        :returns:  An API can be utilized to access a gene list, organism list,
                   taxon-id list, tier list/dict/data-frame, accession
                   list/data-frame, blast query list, mygene information, and
                   missing/duplicate information.
        """

        # Usually only a user would manually add a csv file for their own
        # purposes.
        self.blastn_log.info("Getting the master lists.")
        if csv_file is not None:
            self.__init__(project=self.project, acc_file=csv_file)
            df = self.df
        maf = df
        self.gene_list = maf.index.tolist()
        self.gene_count = len(self.gene_list)

        if self.__hgnc_source is not None:
            self.hgnc_df = self.get_hgnc_gene_info(
                self.gene_list,
                source=self.__hgnc_source,
            )
            self.hgnc_df.to_csv(self.hgnc_path, index=False)

        self.org_list = maf.axes[1].tolist()[1:]
        self.org_count = len(self.org_list)
        self.ncbi_orgs = list(org.replace('_', ' ') for org in self.org_list)

        if self.taxon_file is not None:
            # Load taxon ids from a file
            self.taxon_ids = self.get_file_list(self.taxon_path)
        else:
            # Load taxon ids from a local NCBI taxon database via ete3
            ncbi = NCBITaxa()
            taxon_dict = ncbi.get_name_translator(self.ncbi_orgs)
            self.taxon_ids = list(tid[0] for tid in taxon_dict.values())
            self.taxon_orgs = list(torg for torg in taxon_dict.keys())
            self.taxon_orgs = list(org.replace(' ', '_')
                                   for org in self.taxon_orgs)
            self.taxon_dict = dict(zip(self.taxon_orgs, self.taxon_ids))
            self.taxon_lineage = self.get_taxon_dict()

        self.tier_list = maf['Tier'].tolist()
        self.tier_dict = maf['Tier'].to_dict()
        self.tier_frame_dict = self.get_tier_frame()

        self.acc_dict = self.get_acc_dict()
        self.acc_list = list(self.acc_dict.keys())

        # Get blast query list
        if self.ref_species == 'Homo_sapiens':
            self.blast_human = self.df.Homo_sapiens.tolist()
            self.blast_rhesus = self.df.Macaca_mulatta.tolist()

        # Pre-Blast gene analysis
        if self.__pre_blast is True:
            self.mygene_df = self.blast_utils.my_gene_info(
                acc_dataframe=copy.deepcopy(self.raw_acc_data))
            self.mygene_df.to_csv(self.mygene_path, index=False)

        # Post-Blast accession analysis
        if self.__post_blast:

            # Missing
            self.missing_dict = self.blast_utils.get_miss_acc(
                acc_dataframe=copy.deepcopy(self.raw_acc_data))
            self.missing_genes = self.missing_dict['genes']
            self.missing_gene_count = self.missing_genes['count']
            del self.missing_genes['count']
            self.missing_organsims = self.missing_dict['organisms']
            self.missing_organsims_count = self.missing_organsims['count']
            del self.missing_organsims['count']

            # Duplicates
            duplicate_analysis = self.blast_utils.analyze_duplicate_accessions(
                self.acc_dict,
                self.gene_list,
                self.org_list,
            )
            self.duplicated_dict = duplicate_analysis.groups
            self.duplicated_accessions = self.duplicated_dict['accessions']
            self.duplicated_organisms = self.duplicated_dict['organisms']
            self.duplicated_genes = self.duplicated_dict['genes']
            self.duplicated_random = self.duplicated_dict['random']
            self.duplicated_other = self.duplicated_dict['other']
            self.dup_acc_count = duplicate_analysis.accession_counts
            self.dup_gene_count = duplicate_analysis.gene_counts
            self.dup_org_count = duplicate_analysis.organism_counts

    def get_accession(self, gene, organism):
        """Access a single accession number.

        :param gene:  An input gene.
        :param organism:  An input organism.
        :return:  A single accession number of the target gene/organism.
        """

        maf = self.df
        accession = maf.at[gene, organism]
        if isinstance(accession, float):
            accession = 'missing'
        return accession

    def get_orthologous_gene_sets(self, go_list=None):
        """Access a list of accession numbers.

        :param go_list:  A nested list of gene/organism lists
                         (go_list = [[gene.1, org.1], ... , [gene.n, org.n]]).
                         (Default value = None)
        :return:  An ordered list of accession numbers (or "missing") that
                  correspond to the go_list index.
        """

        if go_list is None:
            accessions = self.acc_list
        else:
            accessions = []
            for gene, organism in go_list:
                accession = self.get_accession(gene, organism)
                accessions.append(accession)
        return accessions

    def get_orthologous_accessions(self, gene):
        """Take a single gene & return a list of accession numbers for the different orthologs.

        :param gene:  An input gene from the accession file.
        :return:  A list of accession numbers that correspond to the orthologs
                  of the target gene.
        """

        maf = self.df
        accession_alignment = maf.T[gene].tolist()[1:]
        return accession_alignment

    def get_tier_frame(self, tiers=None):
        """Organize a dictionary by tier.

        Each tier (key) has a value, which is a data-frame of genes
        associated with that tier.

        :param tiers:  A list of tiers in the accession file.
                       (Default value = None)
        :return:  A nested dictionary for accessing information by tier.
        """

        maf = self.df
        tier_frame_dict = {}
        if tiers is None:
            tiers = maf.groupby('Tier').groups.keys()
        for tier in tiers:
            tier_frame_dict[str(tier)] = maf.groupby('Tier').get_group(tier)
        return tier_frame_dict

    def get_taxon_dict(self):
        """Get the taxonomy information about each organism using ETE3.

        :return:  Returns several dictionaries.  One is a basic organism (key)
        to taxonomy id (value) dictionary, and the other is a lineage
        dictionary with the an organism key and a lineage dictionary as the
        value.  The lineage dictionary keys for each organism are
        ["class", "family", "genus", "kingdowm", "order", "phylum", "species",
        "superkingdom"].
        """
        ncbi = NCBITaxa()
        taxa_dict = {}
        for organism in self.org_list:
            try:
                taxa_dict[organism] = {}
                taxid = self.taxon_dict[organism]
                lineage = ncbi.get_lineage(taxid)
                names = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(lineage)
                for id in lineage:
                    if ranks[id] == 'no rank':
                        continue
                    if ranks[id] not in ['superkingdom', 'kingdom', 'phylum',
                                        'class', 'order', 'family', 'genus', 'species']:
                        continue
                    taxa_dict[organism][ranks[id]] = names[id]
            except KeyError as ke:
                self.blastn_log.exception(ke)
        return taxa_dict

    def get_acc_dict(self):
        """Input a list of accession numbers and return a dictionary with corresponding genes/organisms.

        :return: An accession dictionary who's values are nest gene/organism lists.
        """
        # TODO-ROB set up function to accept a parameter for unique values or
        # potential duplicates

        gene_list = self.gene_list
        org_list = self.org_list
        go = {}
        for gene in gene_list:
            for org in org_list:
                query_acc = self.get_accession(gene, org)
                if query_acc not in go:
                    go[query_acc] = []
                # TODO-ROB: Rework the missing function using this.. maybe??
                elif query_acc == 'missing':
                    continue
                go_list = [gene, org]
                # Append so that duplicates can be identified
                go[query_acc].append(go_list)
        return go


class ComparativeGenetics(BaseComparativeGenetics):
    """Main Comparative Genetics class."""

    def __init__(self, project, template=None, taxon_file=None, ref_species=None,
                 post_blast=False, save_data=True, **kwargs):
        """Inherits BaseComparativeGenetics to build a file layer to the Blast workflow.

        This class handles all of the files before and after the Blast occurs.
        It also uses a building file to start where a previous blast left off.

        :param project:  The name of the project.
        :param template:  A template accession file in the desired format.
                          See the Blast README for an example.
        :param taxon_file:  A list of taxon ids in a text file.
        :param ref_species: A reference species or organism for the blast query.
        :param post_blast:  A flag that triggers the post blast analysis.
        :param save_data:  A flag that indicates whether the data should be
                           saved in an excel file or not.
        :param kwargs:  Mostly used for BaseComparativeGenetics
        :returns:  An API for accessing the various files used before, during,
                   and after blasting."""

        super().__init__(project=project, taxon_file=taxon_file,
                         ref_species=ref_species, post_blast=post_blast,
                         hgnc=False, **kwargs)

        self.postblastlog = LogIt().default(logname="post blast", logfile=None)
        self.ref_species = ref_species

        self.acc_file = template
        # Private variables
        self.__home = os.getcwd()
        if self.taxon_file is not None:
            self.taxon_path = self.project_index / Path(self.taxon_file)
        self.__post_blast = post_blast
        self.save_data = save_data

        if template is not None:
            self.template_filename = template
            self.template_path = self.project_index / \
                Path(self.template_filename)
            self.building_filename = str(template[:-4] + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')
        else:
            self.building_filename = str(self.project + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')

    def add_accession(self, gene, organism, accession):
        """Build an accession file after a Blastn run.

        It finds whether or not the Blast has been interrupted or not, so that
        the Blast can pick up where it left off.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param accession:  The accession of interest.
        :return:
        """

        # TODO-ROB:  Create this in the log file
        if pd.isnull(self.building.at[gene, organism]) is False:
            existing = self.building.at[gene, organism]
            if existing == accession:
                self.blastn_log.warning(self.sep)
                self.blastn_log.warning("Blastn has run on this gene.")
                self.blastn_log.warning("The ACCESSION(%s) for the %s %s gene already exists in our data set."
                                        % (accession, organism, gene))
                self.blastn_log.warning(self.sep)
            else:
                self.blastn_log.critical(self.sep)
                self.blastn_log.critical("BlastN has run on this gene.")
                self.blastn_log.critical("The queried ACCESSION (%s) does not match the existing ACCESSION (%s)"
                                         % (accession, existing))

                self.blastn_log.critical("Queried Accession Gene: %s" % gene)
                self.blastn_log.critical(
                    "Queried Accession Organism: %s" %
                    organism)
                self.blastn_log.critical(
                    "Existing Accession Gene: %s" %
                    self.acc_dict[existing][0][0])
                self.blastn_log.critical(
                    "Existing Accession Organism: %s" %
                    self.acc_dict[existing][0][1])
                self.blastn_log.critical(self.sep)
                self.blastn_log.warning("The existing accession will be overwritten.")

        self.building.at[gene, organism] = accession
        temp = self.building.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        # TODO: make the query organism insert implicit
        temp.insert(2, self.ref_species, self.df[self.ref_species])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_file_path))

    def add_blast_time(
        self,
        gene: str,
        organism: str,
        start: float,
        end: float,
    ) -> None:
        """Build a file that stores the amount of time for each gene to blast.

        This method is similar to the add_accession() method.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param start:  Starting time.
        :param end:  Ending time.
        """
        elapsed_time = end - start
        self.time_dict.setdefault(gene, {})[organism] = elapsed_time
        # Edit the data frame
        self.building_time.at[gene, organism] = elapsed_time
        temp = self.building_time.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        temp.insert(2, self.ref_species, self.df[self.ref_species])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_time_file_path))

    def post_blast_analysis(
        self,
        removed_genes: Sequence[str] | None = None,
    ) -> Path | None:
        """Write duplicate, missing, and removed-gene results to Excel."""
        worksheets: dict[str, pd.DataFrame] = {}

        if removed_genes:
            worksheets["Removed Genes"] = pd.DataFrame(
                {"Removed Genes": list(removed_genes)}
            )

        optional_worksheets = (
            (
                "Duplicate Count by Accession",
                _count_worksheet(self.dup_acc_count),
            ),
            ("Duplicate Count by Gene", _count_worksheet(self.dup_gene_count)),
            (
                "Duplicate Org Groups by Gene",
                _duplicate_group_worksheet(self.duplicated_genes),
            ),
            ("Duplicate Count by Org", _count_worksheet(self.dup_org_count)),
            (
                "Duplicate Gene Groups by Org",
                _duplicate_group_worksheet(self.duplicated_organisms),
            ),
            ("Random Duplicates", _mapping_worksheet(self.duplicated_random)),
            ("Other Duplicates", _mapping_worksheet(self.duplicated_other)),
        )
        worksheets.update(
            {
                sheet_name: worksheet
                for sheet_name, worksheet in optional_worksheets
                if worksheet is not None
            }
        )
        worksheets.update(
            _missing_worksheets(
                self.missing_organsims,
                details_key="missing genes",
                count_sheet_name="Missing Genes Count",
                details_sheet_name="Missing Genes by Org",
            )
        )
        worksheets.update(
            _missing_worksheets(
                self.missing_genes,
                details_key="missing organisms",
                count_sheet_name="Missing Organisms Count",
                details_sheet_name="Missing Organisms by Gene",
            )
        )

        if not worksheets:
            self.postblastlog.warning(
                "Post-BLAST analysis contained no reportable results."
            )
            return None

        output_path = self.data / f"{self.project}_postblastanalysis.xlsx"
        with pd.ExcelWriter(output_path) as workbook:
            for sheet_name, worksheet in worksheets.items():
                worksheet.to_excel(workbook, sheet_name=sheet_name)

        self.postblastlog.info(
            f"Post-BLAST analysis written to {output_path}."
        )
        return output_path
