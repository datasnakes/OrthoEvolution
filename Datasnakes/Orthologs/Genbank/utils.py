from Bio import SeqIO
from pathlib import Path


def multi_fasta_manipulator(full_file, id_file, manipulation='remove', added_name='_G2'):
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
    new_records = set()
    full_file = Path(full_file).open()
    id_file = Path(full_file).open()
    # Turn the id_file into set of ids
    with open(id_file) as id_handle:
        ids = set(record.id for record in SeqIO.parse(id_handle, 'fasta'))
    # Create a new multi-fasta record object using the full_file and the newly generated set of ids
    if manipulation == 'remove':
        new_records.update(set(record for record in SeqIO.parse(full_file, 'fasta') if record.id not in ids))
        print('Sequences have been filtered.')
    # Combine all the FASTA sequence in one record object
    elif manipulation == 'add':
        new_records.update(set(record for record in SeqIO.parse(full_file, 'fasta')))
        new_records.update(set(record for record in SeqIO.parse(id_file, 'fasta')))
        print('Sequences have been added.')

    new_file = Path(full_file).parent / Path(full_file).stem + added_name + Path(full_file).suffix
    SeqIO.write(new_records, new_file, 'fasta')
    print('A new fasta file has been created.')
    return new_file

