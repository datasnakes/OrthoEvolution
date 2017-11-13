# Docs
This folder hosts our tutorial documents used for our readthedocs page as well
as functions to create our docs dynamically.

Our readthedocs page is [here](http://datasnakes-scripts.readthedocs.io/en/master/).

## ❗ Software Dependencies
[Pandoc](http://johnmacfarlane.net/pandoc/) must be installed to use our `PandocConverter` class.

You can install it using the `pypandoc` package as well.

```python
# expects an installed pypandoc: pip install pypandoc
from pypandoc.pandoc_download import download_pandoc
# see the documentation how to customize the installation path
# but be aware that you then need to include it in the `PATH`
download_pandoc()
```

