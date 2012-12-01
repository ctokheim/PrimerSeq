import wx
from wx.lib.pubsub import pub
import threading


class PlotThread(threading.Thread):
    def __init__(self, target, args):
        threading.Thread.__init__(self)
        self.tar = target
        self.args = args
        self.start()

    def run(self):
        output = self.tar(*self.args)  # threaded call
        wx.CallAfter(pub.sendMessage, "plot_update", ())


class UpdateThread(threading.Thread):
    def __init__(self, target, args, update=''):
        threading.Thread.__init__(self)
        self.update = update
        self.tar = target
        self.args = args
        self.start()

    def run(self):
        output = self.tar(*self.args)  # threaded call
        wx.CallAfter(pub.sendMessage, self.update, ())


class RunThread(threading.Thread):
    def __init__(self, target, args, attr='', label='', label_text=''):
        threading.Thread.__init__(self)
        self.label = label
        self.label_text = label_text
        self.tar = target
        self.args = args
        self.attr = attr
        self.start()

    def run(self):
        output = self.tar(*self.args)  # threaded call

        # Only for loading files. Not for when running PrimerSeq.
        if self.attr and self.label and self.label_text:
            wx.CallAfter(pub.sendMessage, "update", ((self.attr, output), (self.label, self.label_text)))
        else:
            wx.CallAfter(pub.sendMessage, "update", (None,))  # need to make this call more elegant

