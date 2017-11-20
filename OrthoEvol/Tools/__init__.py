"""This collection of modules provides useful tools for this package."""


class ToolsWarning(Warning):
    """This is the Tools main warning class.

    The Orthologs package/module will use this warning to alert uses about
    potential tricky code or code under development.

    >>> import warnings
    >>> from OrthoEvol.Tools import ToolsWarning
    >>> warnings.simplefilter('ignore', ToolsWarning)
    """

    pass


class ToolsDevelopmentWarning(ToolsWarning):
    """This is the Tools developmental code warning subclass.

    This warning is for alpha or beta level code which is released as part of
    the standard releases to mark sub-modules or functions for early adopters
    to test & give feedback. Code issuing this warning is likely to change or
    removed in a subsequent release of this package. Such code should NOT be
    used for production/stable code.
    """

    pass
