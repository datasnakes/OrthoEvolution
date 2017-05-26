import webview
import threading
import os
"""
This example demonstrates creating an open file dialog.
"""


def open_file_dialog():
    import time
    time.sleep(5)
    print(
        webview.create_file_dialog(
            webview.OPEN_DIALOG,
            allow_multiple=True,
            directory=os.getcwd()))


if __name__ == '__main__':
    t = threading.Thread(target=open_file_dialog)
    t.start()

    webview.create_window(
        "Open file dialog example",
        "http://www.vallender-genetics.ml/shiny/")
