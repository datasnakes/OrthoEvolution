class OrthoEvolWarning(Warning):
    """This is the OrthoEvol main warning class.

    The OrthoEvol-Scripts package/module will use this warning to alert uses
    about potential tricky code or code under development.

    >>> import warnings
    >>> from OrthoEvol import OrthoEvolWarning
    >>> warnings.simplefilter('ignore', OrthoEvolWarning)
    """

    pass


class OrthoEvolDevelopmentWarning(OrthoEvolWarning):
    """This is the OrthoEvol developmental code warning subclass.

    This warning is for alpha or beta level code which is released as part of
    the standard releases to mark sub-modules or functions for early adopters
    to test & give feedback. Code issuing this warning is likely to change or
    removed in a subsequent release of this package. Such code should NOT be
    used for production/stable code.
    """

    pass


class OrthoEvolDeprecationWarning(OrthoEvolWarning):
    """This is the Deprecation Warning subclass.

    This warning is for code that is no longer maintained and will be removed
    from this project at a later date. It may may not be working as intended,
    and it will not be fixed or edited.
    """

    pass
