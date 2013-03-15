#!/usr/bin/env python
# Copyright (C) 2013 Collin Tokheim
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


class SavePlotsHTML(object):
    """
    SavePlotsHTML creates html for displaying plots about read depth,
    transcript abundance, and primer design. This allows users to view
    the results outside of the PrimerSeq program.
    """
    def __init__(self):
        self.html_content = ''

    def add_heading(self, text, heading_type='h1'):
        """Add a heading element"""
        self.html_content += '<center><%s>%s</%s></center>\n' % (heading_type, text, heading_type)

    def add_line_break(self):
        """Add a br tag"""
        self.html_content += '<br />\n'

    def add_horizontal_line(self):
        """Add a hr tag"""
        self.html_content += '<hr />\n'

    def add_link(self, url, text):
        """Add a hyperlink"""
        self.html_content += '<a href="%s" target="_blank">%s</a>' % (url, text)

    def add_text(self, text):
        """Just add verbatim text to html"""
        self.html_content += text

    def add_img(self, img):
        """Add an image"""
        self.html_content += '<img src="%s" />\n' % img

    def get_html(self):
        """Returns html with all the necessary boilerplate"""
        return """<!DOCTYPE html>
<html>
<head>
    <LINK href="style.css" rel="stylesheet" type="text/css">
    <title>PrimerSeq v1.0.4.beta Output</title>
</head>
<body>
<div class='container'>
%s
</div>
</body>
</html>""" % (self.html_content)

    def __str__(self):
        """Convenience method for get_html"""
        return self.get_html()

    def __repr__(self):
        """Convenience method for printing"""
        return self.get_html()

if __name__ == '__main__':
    """Creates a Hello World html page. Also demonstrates the page style."""
    myHtml = SavePlotsHTML()
    myHtml.add_heading('Hello World!!!')
    print myHtml
