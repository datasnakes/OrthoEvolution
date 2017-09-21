"""Send messages or milestones to a slack channel."""
import configparser
from slacker import Slacker
import os


class Slackify(object):
    def __init__(self, slackconfig='slackconfig.cfg', cfg=True):
        config = configparser.ConfigParser()
        # If there is not a config file
        # Use False if you did not create a config file
        if not cfg:
            apikey = input('Insert your slack apikey here: ')
            if len(apikey) is not 42:  # Standard length of slack apikey
                print('Your slack APIKEY is incorrect.')
                raise ValueError
            self.slack = Slacker(apikey)

        # If there is a config file
        if not os.path.isfile(slackconfig):
            raise OSError
        config.read(slackconfig)
        # HINT Create a config file like the one described in the readme
        apikey = config['APIKEYS']['slack']
        self.slack = Slacker(apikey)

    def upload_img(self, channel, imgfile):
        """Upload an image to a slack channel."""
        self.slack.files.upload(imgfile, channel=channel)

    def upload_file(self, channel, file):
        """Upload files (text/pdf/docx/log) to a slack channel."""
        self.slack.files.upload(file, channel=channel)

    def message_slack(self, channel, message, username):
        """Post a message to slack channel."""
        self.slack.chat.post_message(channel, message, username, as_user=True)