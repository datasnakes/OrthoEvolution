# -*- coding: utf-8 -*-
"""
Date created: Tue Mar  7 17:01:10 2017
Author: S. Hutchins

Script description: Just a few definitions to go with the params_test.py script.

"""

from Bio.Align.Applications import ClustalOmegaCommandline


# Tool that uses the CLustalOmega Command line
def COP1(in_file, out_file, outfmt, logfile):
    """Default parameters with 4 combined iterations"""
    clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, iterations=4, verbose=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP2(in_file, out_file, outfmt, logfile):
    """Default parameters while using a profile with 1 iteration"""
    clustalo_cline = ClustalOmegaCommandline(profile1="profile.fasta", infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, iterations=1, verbose=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP3(in_file, out_file, outfmt, logfile):
    """1 iteration + profile + distmat full iteration"""
    clustalo_cline = ClustalOmegaCommandline(profile1="profile.fasta", infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, iterations=1, verbose=True,
                                             distmat_full_iter=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP4(in_file, out_file, outfmt, logfile):
    """1 iteration + distmat full iteration """
    clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, iterations=1, verbose=True,
                                             distmat_full_iter=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP5(in_file, out_file, outfmt, logfile):
    """Default parameters with 1 interation"""
    clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, iterations=1, verbose=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP6(in_file, out_file, outfmt, logfile):
    """3 combined iterations with 2 hmm iterations"""
    clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, seqtype="DNA", max_hmm_iterations=2,
                                             infmt="fasta", outfmt=outfmt, iterations=3, verbose=True,
                                             threads=8, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;

def COP7(in_file, out_file, outfmt, logfile):
    """Default with auto set to TRUE"""
    clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, seqtype="DNA",
                                             infmt="fasta", outfmt=outfmt, verbose=True,
                                             threads=8, auto=True, force=True, log=logfile)
    stdout, stderr = clustalo_cline()
    clustalo_cline()
    print(stdout, stderr)
    print("\n" + "File has been created." + "\n")
    return;
