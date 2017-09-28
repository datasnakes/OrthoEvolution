slackify
=============
Send updates to Slack about your pipeline's progression!

You can upload a file, image, or send a message to a slack channel once you've
gone through Slack to [generate an API KEY](https://get.slack.help/hc/en-us/articles/215770388-Create-and-regenerate-API-tokens).

The bot you create (if you choose that route) must be invited to the channel you post from the bot in.

After generating an apikey, it's best to create a configuration file so that you can easily keep up with your apikey.
Make sure to practice secure methods. Don't upload your apikey to github as that is very insecure. Keep a local copy of your key.

Example
--------

Message a channel and link to a user with `<@username>` in your message.

```python

slack = Slackify(slackconfig='path/to/slackconfig.cfg')

message_to_channel = 'Hey, <@username>. This is an update for the current script.'

slack.send_msg(channel='channelname', message=message_to_channel)
```