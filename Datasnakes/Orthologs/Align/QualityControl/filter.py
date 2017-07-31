import os
from pathlib import Path
from Datasnakes.Orthologs.Align.QualityControl.guidance2 import Guidance2Commandline
from Datasnakes.Orthologs.Align.QualityControl.pal2nal import Pal2NalCommandline


# For one gene at a time
class QCAlignmentCommandline(Guidance2Commandline, Pal2NalCommandline):

    def __init__(self, na_fasta, aa_fasta, home=os.getcwd(), msaProgram='CLUSTALW', na_bootstraps=10, aa_bootstraps=25,
                 na_seqCutoff=0.6, aa_seqCutoff=0.6, na_colCutoff=0, aa_colCutoff=0.88):
        Guidance2Commandline.__init__(outOrder='as_input', dataset=Path(na_fasta).stem, msaProgram=msaProgram)
        Pal2NalCommandline.__init__(nogap=True, nomismatch=True, output_file=Path(na_fasta).stem + '_P2N.aln')

        # Initialize
        self.home = home
        self.na_fasta = na_fasta
        self.aa_fasta = aa_fasta
        self.gene = Path(self.na_fasta).stem
        self.na_bootstraps = na_bootstraps
        self.aa_bootstraps = aa_bootstraps
        self.na_seqCutoff = na_seqCutoff
        self.aa_seqCutoff = aa_seqCutoff
        self.na_colCutoff = na_colCutoff
        self.aa_colCutoff = aa_colCutoff
        self.msaProgram = msaProgram
        self.na_guidance_path = Path(self.home) / Path('NA_Guidance2')
        self.aa_guidance_path = Path(self.home) / Path('AA_Guidance2')
        Path.mkdir(self.na_guidance_path, exist_ok=True)
        Path.mkdir(self.aa_guidance_path, exist_ok=True)
        # Guidance2 iterations over the nucleic acid sequences
        iteration_flag = True
        while iteration_flag is True:
            self.nucleic_acid_guidance()
        print('')

    def nucleic_acid_guidance(self, seqFile, outDir, program, bootstraps, genCode, seqCutoff, colCutoff, MSA_Param, proc_num):
        seqType = 'nuc'
        print('Align the nucleic acid sequences using Guidance2.')

    @staticmethod
    def amino_acid_guidance(aa_fasta, na_removed_fasta):
        seqType = 'aa'
        print('Align the filtered amino acid sequences using guidance 2')

    def pal2nal_conversion(self, aa_alignment, na_fasta):
        print('Align the nucleic acids using the amino acid alignment.')
