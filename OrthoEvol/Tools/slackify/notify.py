"""Slackify sends messages or milestones to a slack channel."""
import configparser
import os
from slacker import Slacker
import slacker


class Slackify(object):
    """Send messages or milestones to a slack channel."""
    def __init__(self, slackconfig='slackconfig.cfg', cfg=True):
        """Configure Slackify.

        :param slackconfig: Path to config file.
        :param cfg: Flag that is default True.
        """
        config = configparser.ConfigParser()
        # If there is not a config file
        # Use False if you did not create a config file
        if not cfg:
            apikey = input('Insert your slack apikey here: ')
            if len(apikey) is not 42:  # Standard length of slack apikey
                raise ValueError('Your slack APIKEY is incorrect.')

            slack = Slacker(apikey)

        else:
            # If there is a config file
            if not os.path.isfile(slackconfig):
                raise FileNotFoundError('Slack configuriation file not found.')
            config.read(slackconfig)
            # HINT Create a config file like the one described in the readme
            apikey = config['APIKEYS']['slack']
            slack = Slacker(apikey)

        self.slack = slack

    def _get_channel_id(self, channel):
        """Get a channel id for uploading files.

        :param channel: Name of the channel to get an id for.
        """

        channel_id = self.slack.channels.get_channel_id(channel)
        return channel_id

    def upload_file(self, file, channel):
        """Upload files (text/pdf/docx/log/image) to a slack channel.

        :param file: Path to a file to upload.
        :param channel: Name of the channel to upload the file to.
        """

        channel_id = self._get_channel_id(channel)
        self.slack.files.upload(file_=file, channels=channel_id)

    def send_msg(self, channel, message):
        """Post a message to slack channel.

        Send a message to a user using <@username>

        :param channel: Name of the channel to send a message to.
        :param message: Message to send to a channel.
        """

        # With as user as True, the predefined bot name is used
        try:
            self.slack.chat.post_message(channel, message, as_user=True)
            print('Your message was sent.')
        except slacker.Error:
            print('Your message was not sent!.')
            # TODO add a traceback here.
            # TODO also use contextlib

    def list_users(self):
        """List all users for your slack organization."""

        response = self.slack.users.list()
        users = [username['name'] for username in response.body['members']]
        return users

    def list_channels(self):
        """List all channels for your slack organization."""

        response = self.slack.channels.list()
        channels = [channel['name'] for channel in response.body['channels']]
        return channels

    def log2slack(self):
        """Send a formatted text string to slack similar to logging."""

        raise NotImplementedError('This function is not yet implemented.')
        # TODO One day...create logging format for logging to slack.
