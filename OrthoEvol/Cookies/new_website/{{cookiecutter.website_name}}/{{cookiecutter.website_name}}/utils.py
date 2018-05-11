# -*- coding: utf-8 -*-
"""Helper utilities and decorators."""
from flask import flash


def flash_errors(form, category='warning'):
    """Flash all errors for a form.

    :param form: 
    :param category:  (Default value = 'warning')

    """
    for field, errors in form.errors.items():
        for error in errors:
            flash('{0} - {1}'.format(getattr(form, field).label.text, error), category)
