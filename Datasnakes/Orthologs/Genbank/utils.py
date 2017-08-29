from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
import os
from pathlib import Path
from tempfile import TemporaryFile
import itertools


def multi_fasta_manipulator(target_file, reference, output, manipulation='remove'):
    # Inspired by the BioPython Tutorial and Cookbook ("20.1.1 Filtering a sequence file")
    """
    This method manipulated selected sequences in a multi-FASTA files.  The original
    purpose was to filter files created by the GUIDANCE2 alignment program, but
    the function has been made in order to accommodate other situations as well.

    :param target_file:  Target multi-FASTA file.
    :param reference:  Selected sequences for removal in a multi-FASTA file.
    :param manipulation:  Type of manipulation.  (remove, add, tbd..)
    :param added_name:  The output file uses this parameter to name itself.
    :return:  A multi-FASTA file with filter sequences.
    """
    # Create path variables
    file_name = output
    new_file = Path(target_file).parent / Path(file_name)
    # Create a new multi-fasta record object using the target_file, reference, and output
    # Remove specific sequences from a fasta file
    if manipulation is 'remove':
        multi_fasta_remove(target_file, reference, new_file)
    # Combine all the FASTA sequence in one record object
    elif manipulation is 'add':
        muli_fasta_add(target_file, reference, new_file)
    # Sort one fasta file based on the order of another
    # Works for alignments
    elif manipulation is 'sort':
        multi_fasta_sort(target_file, reference, new_file)

    print('A new fasta file has been created.')
    return new_file


def dir_config(path, tier_frame_dict):
    """
    Configure the genbank directories.
    :param path: Path to create directory structure.
    :param tier_frame_dict:  Dictionary from the blast super class.
    :return:  Creates a directory structure as follows
        --Tier_1
            --Gene_1
            --Gene_M
        --Tier_N
            --Gene_M+1
            --Gene_N
    """
    for G_KEY in tier_frame_dict.keys():
        tier = G_KEY
        tier_path = path / Path(tier)
        Path.mkdir(tier_path, parents=True, exist_ok=True)
        for GENE in tier_frame_dict[tier].T:
            gene_path = tier_path / Path(GENE)
            Path.mkdir(gene_path)


def multi_fasta_remove(target_file, reference, new_file):
    rem_file = new_file.stem + '_removed' + new_file.suffix
    rem_file = new_file.parent / Path(rem_file)
    # Turn the reference_file into set of ids
    if os.path.isfile(reference):
        ids = set(record.id for record in SeqIO.parse(reference, 'fasta'))
    elif isinstance(reference, list):
        ids = reference

    new_records = (record for record in SeqIO.parse(target_file, 'fasta') if record.id not in ids)
    old_records = (record for record in SeqIO.parse(target_file, 'fasta') if record.id in ids)

    print('Sequences have been filtered.')
    SeqIO.write(new_records, str(new_file), 'fasta')
    SeqIO.write(old_records, str(rem_file), 'fasta')


def muli_fasta_add(target_file, reference, new_file):
    # TODO-ROB:  Check for duplicates.
    # Concatenate the multifasta files together by chaining the SeqIO.parse generators
    # Allows one to overwrite a file by using temporary files for storage
    # adding generators - https://stackoverflow.com/questions/3211041/how-to-join-two-generators-in-python
    if os.path.isfile(reference):
        with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
            new_records = itertools.chain(SeqIO.parse(target_file, 'fasta', ), SeqIO.parse(reference, 'fasta'))
            count = SeqIO.write(new_records, tmp_file, 'fasta')
            tmp_file.seek(0)
            print('temp file count: ' + str(count))
            SeqIO.write(SeqIO.parse(tmp_file, 'fasta'), str(new_file), 'fasta')
        print('Sequences have been added.')
    else:
        print('You can only add files together.  Not python objects.')


def multi_fasta_sort(target_file, reference, new_file):
    # TODO-ROB:  Check for duplicates.
    with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
        aln = MultipleSeqAlignment([])

    # For a reference fasta file make a tuple from ids
    if os.path.isfile(reference):
        sorted_list = []
        for s in SeqIO.parse(reference, 'fasta'):
            sorted_list.append(s.id)
        sorted_handle = tuple(sorted_list)
    # For a reference list port in as a tuple
    elif isinstance(reference, list):
        sorted_handle = tuple(reference)

    # Parese the reference tuple above the unsorted file
    for sorted_seq_record in sorted_handle:
        for unsorted_aln_record in SeqIO.parse(target_file, 'fasta'):
            # If an appropriate id is found, then append to the MSA object.
            if unsorted_aln_record.id == sorted_seq_record.id:
                print(unsorted_aln_record.id)
                print(sorted_seq_record.id)
                aln.append(unsorted_aln_record)  # MSA object
                break
    count = AlignIO.write(aln, tmp_file, 'fasta')
    tmp_file.seek(0)
    print('temp file count: ' + str(count))
    AlignIO.write(AlignIO.read(tmp_file, 'fasta'), str(new_file), 'fasta')
    print('Alignment has been sorted.')

