# def git_ignore(self, path):
#     """Get the ignored file patterns from the .gitignore file in the repo."""
#     with open(Path(path) / Path('.gitignore'), 'r', newline='') as ignore:
#         ignored = ignore.read().splitlines()
#     return ignored
#
# # Map the main project directory.
# def get_dir_map(self, top, gitignore=None):
#     # TODO-ROB:  Change ignore to a .gitignore filename
#     default_ignore = self.git_ignore(top)
#     if gitignore is not None:
#         gitignore += default_ignore
#     else:
#         gitignore = default_ignore
#     # Treelib will help to map everything and create a json at the same time
#     tree = Tree()
#     tree.create_node('.', top)
#     # Walk through the directory of choice (top)
#     # Topdown is true so that directories can be modified in place
#     for root, dirs, files in os.walk(top, topdown=True):
#         # Only remove directories from the top
#         if root == top:
#             print(root)
#             try:
#                 dirs[:] = set(dirs) - set(gitignore)  # Remove directories from os.walk()
#                 print(dirs)
#             except ValueError:
#                 pass
#         for d in dirs:
#             rd = str(Path(root) / Path(d))
#             tree.create_node(d, identifier=rd, parent=root)
#         for f in files:
#             tree.create_node(f, parent=root)
#     return tree

# def get_newick_dir_map(self, top, ignore=None):
#     """Takes a treelib tree created by get_dir_map and returns
#     a tree a dir_map in newick format.  This will be useful for Bio.Phylo
#     applications."""
#     """
#     :param top (path):  The root at which the directory map is made.
#     :param ignore (list):  The files to ignore.  The  get_dir_map function
#     adds this to the .gitignore list.
#     :return (tree):  A newick formatted string in style #8.  Can be used
#     with the ete3.Tree() class.
#     """

    # tree = Tree()
    # t = tree.get_dir_map(top, ignore)
    # Ntree = tree.parse_newick_json()
    # return Ntree

# def get_ete3_tree(self, top, tree=None):
#     """Create the ete3 tree."""
#     if not tree:
#         tree = self.get_newick_dir_map(top)
#     t = ete3.Tree(tree, format=8)
#     return t

    # DEPRECATED Change this IN OTHER CLASSES
    # def path_list_make(self, path, o_path=None):
    # Takes a path and reduces it to a list of directories within the project
    # An optional attribute (o_path) is give so that a deeper path within
    # the project can be used

    # //TODO-ROB utilize Path.mkdir(parents=TRUE) instead
    # DEPRECATED Change this IN OTHER CLASSES
    #  def dir_make(self, path, path_list):
    # Takes a path list which is a list of folder names
    # path_list created by str(path).split('/')
    # The path_list appends to path, which is already an established
    # directory inside the project

    # # //TODO-ROB Change to using a compression module https://pymotw.com/2/compression.html
    # DEPRECATED Change this IN OTHER CLASSES
    # def dir_archive(self, path, path_list):
    #     # Use the path that you want to update/add to
    #     # Returns path and the time stamp (could be None)