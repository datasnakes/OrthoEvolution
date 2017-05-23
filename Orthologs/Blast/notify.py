from OrthoTools import message_slack
import argparse
import textwrap
import logging as log
from datetime import datetime as d

#------------------------------------------------------------------------------
# Set up the logger
format1 = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
format2 = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives

# Set up the blastn logger & log file
LOGFORMAT = '%(name)s - [%(levelname)-2s]: %(message)s'
log.basicConfig(level=log.DEBUG,
                format=LOGFORMAT,
                filename="logs/slack_notify.log")
slack_log = log.getLogger('Slack Message')

#------------------------------------------------------------------------------
__author__ = 'SDH'

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
                                    This is a command line wrapper for the Vallender Lab's Slack channel.

                                    Channels = #hall-project, #drd4-project, #karg-project,
                                               #l' '''))
parser.add_argument("-c", "--channel", help="Input a channel name", required=True)
parser.add_argument("-m", "--message", help="Write a message here", required=True)
parser.add_argument("-u", "--username", help="Input a username", required=True)
args = parser.parse_args()



message_slack(args.channel, args.message, args.username)
print('Your message was posted to Slack.')
slack_log.info('You posted to the %s channel with the user, %s.' % (args.channel, args.username))
log.shutdown()
