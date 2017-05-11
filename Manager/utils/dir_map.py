import treelib
from pathlib import Path
import os


def git_ignore(path):
    """Get the ignored file patterns from the .gitignore file in the repo."""
    with open(Path(path) / Path('.gitignore'), 'r', newline='') as ignore:
        ignored = ignore.read().splitlines()
    return ignored


def create_newick_node(dir_map, length, tag=None, identifier=None, parent=None, data=None):
    bl = {}
    bl['name'] = tag
    # length must be an int
    bl['branch_length'] = int(length)
    dir_map.create_node(tag, identifier, parent, data=bl)



def to_jsonnewick(dir_map, nid=None, key=None, sort=True, reverse=False, with_data=False):
    # For True newick files add a data parameter for branch length
    # ETE3 format #7 newick tree file format
    if nid is None:
        nid = dir_map.root # identifier (path)
    ntag = dir_map[nid].tag # folder name
    tree_dict = {
        'name': ntag,
        'children': []}
    if with_data is True:
        tree_dict['children'].append(dir_map[nid].data)
    if dir_map[nid].expanded:
        queue = [dir_map[i] for i in dir_map[nid].fpointer]
        if key is None:
            key = (lambda x: x)
        if sort:
            queue.sort(key=key, reverse=reverse)
        for elem in queue:
            tree_dict['children'].append(
                to_jsonnewick(elem.identifier, with_data=with_data, sort=sort, reverse=reverse))
        if len(tree_dict['children']) == 0:
            if not with_data:
                tree_dict = {"name": dir_map[nid].tag}
            else:
                tree_dict = {"branch_length": dir_map[nid].data}
        return tree_dict


def _parse_json(json_obj):
    """Edited from Source:  https://github.com/okeeffdp/json_to_newick"""
    try:
        # Test is the key 'name' in the current level of the dictionary.
        newick = json_obj['name']
    except KeyError:
        # Catch no 'name' trees and start the newick string with empty qoutes
        newick = ''

    if 'branch_length' in json_obj:
        newick = newick + ':' + str(json_obj['branch_length'])
    # If there are 'children'
    if 'children' in json_obj:
        # Initialise a list to contain the daughter info
        info = []
        # For each child, treat it as a new dictionary object
        for child in json_obj['children']:
            # parse the newick string straight into it the list
            info.append(_parse_json(child))
        # join all the daughter info together with a comma
        info = ','.join(info)
        # Concatenate all the children together at the start of the parent newick string
        newick = '(' + info + ')' + newick
    return newick


# Map the main project directory.
def get_dir_map(top, gitignore=None):
    # TODO-ROB:  Change ignore to a .gitignore filename
    default_ignore = git_ignore(os.getcwd())
    if gitignore is not None:
        gitignore += default_ignore
    else:
        gitignore = default_ignore
    # Treelib will help to map everything and create a json at the same time
    tree = treelib.Tree()
    tree.create_node('.', top)
    # Walk through the directory of choice (top)
    # Topdown is true so that directories can be modified in place
    for root, dirs, files in os.walk(top, topdown=True):
        # Only remove directories from the top
        if root == top:
            print(root)
            try:
                dirs[:] = set(dirs) - set(gitignore)  # Remove directories from os.walk()
                print(dirs)
            except ValueError:
                pass
        for d in dirs:
            rd = str(Path(root) / Path(d))
            tree.create_node(d, identifier=rd, parent=root)
        for f in files:
            tree.create_node(f, parent=root)
    return tree


def get_newick_dir_map(top, ignore=None):
    """Takes a treelib tree created by get_dir_map and returns
    a tree a dir_map in newick format.  This will be useful for Bio.Phylo
    applications."""

    tree = get_dir_map(top, ignore)
    Ntree = _parse_json(to_jsonnewick(tree))
    return Ntree