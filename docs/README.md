# Docs

This folder hosts our tutorial documents used for our readthedocs page as well
as functions to create our docs dynamically.

Our readthedocs page is [here](http://orthoevolution.readthedocs.io/en/master/).

## Software Dependencies

[Pandoc](http://johnmacfarlane.net/pandoc/) must be installed to use our `PandocConverter` class.

You can install it using the `pypandoc` package as well.


`pip install pypandoc`

```python
# expects an installed pypandoc: pip install pypandoc
from pypandoc.pandoc_download import download_pandoc
# see the documentation how to customize the installation path
# but be aware that you then need to include it in the `PATH`
download_pandoc()
```

## How-To: Create Our Docs
1. Run createdocs.py using `python createdocs.py` to regenerate docs.
2. Edit docs for specific add ins.
For example, for main documentation in subfolders such as orthologs, reference
the documentation files of submodules using `submodule <submodulereadme.rst>`__
3. Test the docs by running `make html` in the `source` folder.
4. Commit the changes and push to the branch.

### Creating the modules directory for apidocs
Perform the below command in the root directory of this package. First, remove
all of the existing files.
```bash
rm -rf Docs/docs/source/modules/*.rst

sphinx-apidoc OrthoEvol/ -o Docs/docs/source/modules
```

