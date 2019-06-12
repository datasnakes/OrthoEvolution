from .qsub import BaseQsub, Qsub

__all__ = ("BaseQsub",
           "Qsub",
           "_setup_yaml",
           "NotRequiredIf",
           )

def represent_dictionary_order(cls, dict_data):
    return cls.represent_mapping('tag:yaml.org,2002:map', dict_data.items())


def _setup_yaml():
    """ https://stackoverflow.com/a/8661021 """
    yaml.add_representer(OrderedDict, represent_dictionary_order)


# Borrowed from the discussion in the link below:
# https://stackoverflow.com/questions/44247099/click-command-line-interfaces-make-options-required-if-other-optional-option-is
class NotRequiredIf(click.Option):
    def __init__(self, *args, **kwargs):
        self.not_required_if = kwargs.pop('not_required_if')
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs['help'] = (kwargs.get('help', '') +
                          ' NOTE: This argument is mutually exclusive with %s.' %
                          self.not_required_if
                          ).strip()
        super(NotRequiredIf, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        we_are_present = self.name in opts
        other_present = self.not_required_if in opts

        if other_present:
            if we_are_present:
                raise click.UsageError(
                    "Illegal usage: `%s` is mutually exclusive with `%s`." % (
                        self.name, self.not_required_if))
            else:
                self.prompt = None

        return super(NotRequiredIf, self).handle_parse_result(
            ctx, opts, args)