from Bio import SeqIO
from pathlib import Path
from tempfile import TemporaryFile
import itertools

def multi_fasta_manipulator(full_file, id_file, output, manipulation='remove'):
    # Inspired by the BioPython Tutorial and Cookbook ("20.1.1 Filtering a sequence file")
    """
    This method manipulated selected sequences in a multi-FASTA files.  The original
    purpose was to filter files created by the GUIDANCE2 alignment program, but
    the function has been made in order to accommodate other situations as well.

    :param full_file:  Target multi-FASTA file.
    :param id_file:  Selected sequences for removal in a multi-FASTA file.
    :param manipulation:  Type of manipulation.  (remove, add, tbd..)
    :param added_name:  The output file uses this parameter to name itself.
    :return:  A multi-FASTA file with filter sequences.
    """
    # Create path variables
    file_name = output
    new_file = Path(full_file).parent / Path(file_name)
    # Turn the id_file into set of ids
    ids = set(record.id for record in SeqIO.parse(id_file, 'fasta'))
    # Create a new multi-fasta record object using the full_file and the newly generated set of ids
    if manipulation is 'remove':
        new_records = (record for record in SeqIO.parse(full_file, 'fasta') if record.id not in ids)
        print('Sequences have been filtered.')
        SeqIO.write(new_records, str(new_file), 'fasta')
    # Combine all the FASTA sequence in one record object
    elif manipulation == 'add':
        # Concatenate the multifasta files together by chaining the SeqIO.parse generators
        # Allows one to overwrite a file by using temporary files for storage
        # adding generators - https://stackoverflow.com/questions/3211041/how-to-join-two-generators-in-python
        with TemporaryFile('r+', dir=Path(full_file).parent) as tmp_file:
            new_records = itertools.chain(SeqIO.parse(full_file, 'fasta',), SeqIO.parse(id_file, 'fasta'))
            count = SeqIO.write(new_records, tmp_file, 'fasta')
            tmp_file.seek(0)
            print('temp file count: ' + str(count))
            SeqIO.write(SeqIO.parse(tmp_file, 'fasta'), str(new_file), 'fasta')
        print('Sequences have been added.')
    print('A new fasta file has been created.')
    return new_file

