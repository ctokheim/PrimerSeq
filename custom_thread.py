#!/usr/bin/env python
# Copyright (C) 2012  Collin Tokheim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import wx
from wx.lib.pubsub import pub
import threading

# handle uncaught exception imports
import sys


class PlotThread(threading.Thread):
    def __init__(self, target, args):
        threading.Thread.__init__(self)
        self.tar = target
        self.args = args
        self.start()

    def run(self):
        try:
            output = self.tar(*self.args)  # threaded call
            wx.CallAfter(pub.sendMessage, "plot_update", ())
        except:
            wx.CallAfter(pub.sendMessage, "plot_error", ())


class UpdateThread(threading.Thread):
    def __init__(self, target, args, update=''):
        threading.Thread.__init__(self)
        self.update = update
        self.tar = target
        self.args = args
        self.start()

    def run(self):
        try:
            output = self.tar(*self.args)  # threaded call
            wx.CallAfter(pub.sendMessage, self.update, ())
        except:
            wx.CallAfter(pub.sendMessage, "update_after_error", (None,))


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
        try:
            output = self.tar(*self.args)  # threaded call

            # Only for loading files. Not for when running PrimerSeq.
            if self.attr and self.label and self.label_text:
                wx.CallAfter(pub.sendMessage, "update", ((self.attr, output), (self.label, self.label_text)))
            else:
                wx.CallAfter(pub.sendMessage, "update", (None,))  # need to make this call more elegant
        except:
            wx.CallAfter(pub.sendMessage, "update_after_error", (None,))  # need to make this call more elegant


class RunPrimerSeqThread(threading.Thread):
    def __init__(self, target, args, attr='', label='', label_text=''):
        threading.Thread.__init__(self)
        self.label = label
        self.label_text = label_text
        self.tar = target
        self.args = args
        self.attr = attr
        self.start()

    def run(self):
        try:
            output = self.tar(*self.args)  # threaded call
            wx.CallAfter(pub.sendMessage, "update_after_run", (self.args[0]['output'],))  # need to make this call more elegant
        except:
            wx.CallAfter(pub.sendMessage, "update_after_error", (None,))  # need to make this call more elegant


