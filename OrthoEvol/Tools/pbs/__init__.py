from .qsub import BaseQsub, Qsub
from .qstat import BaseQstat, Qstat, MultiQstat, TargetJobKeyError


__all__ = ("BaseQsub",
           "Qsub",
           "BaseQstat",
           "Qstat",
           "MultiQstat",
           "TargetJobKeyError"
           )


# # Borrowed from the discussion in the link below:
# # https://stackoverflow.com/questions/44247099/click-command-line-interfaces-make-options-required-if-other-optional-option-is
# class NotRequiredIf(click.Option):
#     def __init__(self, *args, **kwargs):
#         self.not_required_if = kwargs.pop('not_required_if')
#         assert self.not_required_if, "'not_required_if' parameter required"
#         kwargs['help'] = (kwargs.get('help', '') +
#                           ' NOTE: This argument is mutually exclusive with %s.' %
#                           self.not_required_if
#                           ).strip()
#         super(NotRequiredIf, self).__init__(*args, **kwargs)
#
#     def handle_parse_result(self, ctx, opts, args):
#         we_are_present = self.name in opts
#         other_present = self.not_required_if in opts
#
#         if other_present:
#             if we_are_present:
#                 raise click.UsageError(
#                     "Illegal usage: `%s` is mutually exclusive with `%s`." % (
#                         self.name, self.not_required_if))
#             else:
#                 self.prompt = None
#
#         return super(NotRequiredIf, self).handle_parse_result(
#             ctx, opts, args)
