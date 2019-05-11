"""Optimized for use with local/standalone NCBI BLAST 2.8.1."""
from OrthoEvol.Orthologs.Blast.base_blastn import BaseBlastN


class OrthoBlastN(BaseBlastN):
    """Combines Project Management features with NCBI's Blast+."""

    def __init__(self, project="orthology-inference", method=3, template=None,
                 save_data=True, **kwargs):
        """This class inherits from the CompGenFiles class.

        This class utilizes it's parent classes to search a standalone
        Blast database for specific orthologs of a gene using a query organism
        (usually human).  The best hits from the Blast are filtered for the
        best option in order to get the most accuarate accession numbers for
        downstream analysis.

        :param project:  The project name (Default: 'orthologs')
        :param method: Method used for blasting. (Default: 3)
        :param template:  The accession file template.
        :param save_data:  A flag for saving the post_blast data to an excel file.
        :param kwargs:
        """
        super().__init__(project=project, method=method, template=template, save_data=save_data, **kwargs)

    def run(self):
        self.configure(self.blast_human, 'Homo_sapiens', auto_start=True)
