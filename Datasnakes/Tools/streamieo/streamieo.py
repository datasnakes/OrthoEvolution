from subprocess import Popen, PIPE
from threading import Thread
from queue import Queue, Empty
io_q = Queue()

# Implemented from http://www.sharats.me/posts/the-ever-useful-and-neat-subprocess-module/
# In the "Watching both stdout and stderr" section


def stream_watcher(identifier, stream):

    for line in stream:
        io_q.put((identifier, line))

    if not stream.closed:
        stream.close()


proc = Popen('svn co svn+ssh://myrepo', stdout=PIPE, stderr=PIPE)

Thread(target=stream_watcher, name='stdout-watcher',
       args=('STDOUT', proc.stdout)).start()
Thread(target=stream_watcher, name='stderr-watcher',
       args=('STDERR', proc.stderr)).start()


def printer():
    while True:
        try:
            # Block for 1 second.
            item = io_q.get(True, 1)
        except Empty:
            # No output in either streams for a second. Are we done?
            if proc.poll() is not None:
                break
        else:
            identifier, line = item
            print(identifier + ':', line)


Thread(target=printer, name='printer').start()