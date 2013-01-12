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
        self.css_rules = ''

    class CssProperty(object):
        """CssProperty handles construction of a single id/class css string"""
        def __init__(self, name, properties):
            self.name = name  # .myClassName
            self.properties = properties  # [['margin', 0], ...]
            self.css = ''  # css contents to return
            self.set_css_rules()

        def set_css_rules(self):
            """Creates css for a single id/class"""
            self.css = '%s{' % self.name
            for prop_name, prop_value in self.properties:
                self.css.append('%s:%s;' % (str(prop_name), str(prop_value)))
            self.css.append('}\n')

        def get_css(self):
            """Css string getter"""
            return self.css

        def __str__(self):
            self.get_css()

    def set_css(self):
        """Construct the necessary css rules"""
        self.css_rules += str(self.CssProperty('body, html',
            [['margin', 0],
            ['padding', 0]]))
        self.css_rules += str(self.CssProperty('.container',
            [['margin-left', 'auto'],
            ['margin-right', 'auto'],
            ['width', '800px']]))

    def get_html(self):
        return """<html>
<head>
    <title>PrimerSeq v1.0.0 Output</title>
</head>
<body>
<div class='container'>
%s
</div>
</body>
%s
</html>""" % (self.html_content, self.css_rules)
