"""This collection of modules is ideal for dealing with orthologous genes in
Python. It has been built upon BioPython and ETE3 and provides a great blueprint
for utilizing BioPython and ETE3 to achieve analysis of biological data.
"""


class OrthologsWarning(Warning):
    """This is the Orthologs main warning class.

    The Orthologs package/module will use this warning to alert uses about
    potential tricky code or code under development.

    >>> import warnings
    >>> from datasnakes.Orthologs import OrthologsWarning
    >>> warnings.simplefilter('ignore', OrthologsWarning)
    """

    pass


class OrthologsDevelopmentWarning(OrthologsWarning):
    """This is the Orthologs developmental code warning subclass.

    This warning is for alpha or beta level code which is released as part of
    the standard releases to mark sub-modules or functions for early adopters
    to test & give feedback. Code issuing this warning is likely to change or
    removed in a subsequent release of this package. Such code should NOT be
    used for production/stable code.
    """

    pass
