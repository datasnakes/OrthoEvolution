# PhyloTree Documentation

PhlyoTree is a simple and useful module to help quickly view and create phylogenetic
trees from existing tree files.

## Example

### Draw and save a newick formatted tree

```python
from OrthoEvol.Orthologs.Phylogenetics.TreeViz import TreeViz

t = TreeViz(path='tree.txt', tree_format='newick')

t.draw_tree()
t.save_tree('example.png')
```

## Notes

THIS MODULE IS UNDER DEVELOPMENT!!!!