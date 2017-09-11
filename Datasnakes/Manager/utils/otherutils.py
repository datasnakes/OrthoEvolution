"""Other utilities optimized for the Orthologs Project."""
import pandas as pd


def SplitList(listname, basefilename, n):
    """Split a long list into chunks and save chunks as a text file."""
    # Split the list into chunks
    chunks = [listname[x:x + n] for x in range(0, len(listname), n)]
    list_group = []
    num_lists = len(chunks)

    # Name and save the lists
    for chunk, num in zip(chunks, range(0, num_lists)):
        listdf = pd.DataFrame(chunk)
        n = basefilename + '_list_' + str(num)
        listdf.to_csv(n + ".txt", index=False, header=None)
        list_group.append(n)
    return list_group


def FormatList(input_list):
    """Remove spaces from list items and turn those spaces into underscores."""
    output_list = []
    for item in input_list:
        item = str(item)
        item = item.replace(" ", "_")
        output_list.append(item)
        return output_list
