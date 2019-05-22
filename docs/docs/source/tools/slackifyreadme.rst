slackify
========

Send updates to Slack about your pipeline's progression!

You can upload a file, image, or send a message to a slack channel once
you've gone through Slack to `generate an API
KEY <https://get.slack.help/hc/en-us/articles/215770388-Create-and-regenerate-API-tokens>`__.

The bot you create (if you choose that route) must be invited to the
channel you post from the bot in.

After generating an apikey, it's best to create a configuration file so
that you can easily keep up with your apikey. Make sure to practice
secure methods. Don't upload your apikey to github as that is very
insecure. Keep a local copy of your key.

Examples
--------

Import the class and set up the slack handler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.slackify import Slackify

    slack = Slackify(slackconfig='path/to/slackconfig.cfg')

Your config file should look as such:

.. code:: python

    [APIKEYS]
    slack = apikeystring

Message a channel and link to a user with ``<@username>`` in your message
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    message_to_channel = 'Hey, <@username>. This is an update for the current script.'

    slack.send_msg(channel='channelname', message=message_to_channel)

Get all users and channels
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    slack.list_users() # Returns a list of all users.
    slack.list_channels() # Returns a list of channels

Upload a file
~~~~~~~~~~~~~

The file can be an image, pdf, doc, text, python file, etc

.. code:: python

    slack.upload_file()
