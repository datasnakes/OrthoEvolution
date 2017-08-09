from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from pathlib import Path
from tempfile import TemporaryFile
import itertools

def multi_fasta_manipulator(target_file, reference_file, output, manipulation='remove'):
    # Inspired by the BioPython Tutorial and Cookbook ("20.1.1 Filtering a sequence file")
    """
    This method manipulated selected sequences in a multi-FASTA files.  The original
    purpose was to filter files created by the GUIDANCE2 alignment program, but
    the function has been made in order to accommodate other situations as well.

    :param target_file:  Target multi-FASTA file.
    :param reference_file:  Selected sequences for removal in a multi-FASTA file.
    :param manipulation:  Type of manipulation.  (remove, add, tbd..)
    :param added_name:  The output file uses this parameter to name itself.
    :return:  A multi-FASTA file with filter sequences.
    """
    # Create path variables
    file_name = output
    new_file = Path(target_file).parent / Path(file_name)
    if reference_file is not None:
        # Turn the reference_file into set of ids
        ids = set(record.id for record in SeqIO.parse(reference_file, 'fasta'))
    # Create a new multi-fasta record object using the target_file and the newly generated set of ids
    if manipulation is 'remove':
        new_records = (record for record in SeqIO.parse(target_file, 'fasta') if record.id not in ids)
        print('Sequences have been filtered.')
        SeqIO.write(new_records, str(new_file), 'fasta')
    # Combine all the FASTA sequence in one record object
    elif manipulation is 'add':
        # Concatenate the multifasta files together by chaining the SeqIO.parse generators
        # Allows one to overwrite a file by using temporary files for storage
        # adding generators - https://stackoverflow.com/questions/3211041/how-to-join-two-generators-in-python
        with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
            new_records = itertools.chain(SeqIO.parse(target_file, 'fasta', ), SeqIO.parse(reference_file, 'fasta'))
            count = SeqIO.write(new_records, tmp_file, 'fasta')
            tmp_file.seek(0)
            print('temp file count: ' + str(count))
            SeqIO.write(SeqIO.parse(tmp_file, 'fasta'), str(new_file), 'fasta')
        print('Sequences have been added.')
    elif manipulation is 'sort':
        with TemporaryFile('r+', dir=str(Path(target_file).parent)) as tmp_file:
            aln = MultipleSeqAlignment([])
            for sorted_seq_record in SeqIO.parse(reference_file, 'fasta'):
                for unsorted_aln_record in SeqIO.parse(target_file, 'fasta'):
                    if unsorted_aln_record.id == sorted_seq_record.id:
                        print(unsorted_aln_record.id)
                        print(sorted_seq_record.id)
                        aln.append(unsorted_aln_record)
                        break
            count = AlignIO.write(aln, tmp_file, 'fasta')
            tmp_file.seek(0)
            print('temp file count: ' + str(count))
            AlignIO.write(AlignIO.read(tmp_file, 'fasta'), str(new_file), 'fasta')
        print('Alignment has been sorted.')
    print('A new fasta file has been created.')
    return new_file

