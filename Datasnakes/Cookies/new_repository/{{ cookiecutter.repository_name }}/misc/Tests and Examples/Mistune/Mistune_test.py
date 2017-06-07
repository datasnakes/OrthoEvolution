# -*- coding: utf-8 -*-
"""
@author: S. Hutchins

Mistune is a markdown parser in pure Python with renderer features.

"""
# Import the modules
import mistune

# Create a variable for mistune.markdown
m = mistune.markdown

# Create a markdown file
l1 = m('This is a test for **MISTUNE**!')
# output: <p>I am using <strong>markdown</strong></p>

l2 = m('*You are welcome!* :sparkles:')
# output: <p><em>You are welcome!</em> :sparkles:</p>

l3 = m('Read the Docs: [Mistune](http://mistune.readthedocs.io/en/latest/)')
# output: <p>Read the Docs: <a href="http://mistune.readthedocs.io/en/latest/">Mistune</a></p>

with open('ReadMe.md', 'w') as mdfile:
    mdfile.writelines(l1 + "\n" + l2 + "\n" + l3)