'''
Created on Feb 8, 2012
Intended to significantly aid in the drawing of exon structure and text box positioning.
Shape objects are intended to store a lot of attributes about a shape so that changing a shapes
position becomes easy when calling a method. The JunctionLine class acts to create lines and text
boxes necessary to effectively show the Alternative Splicing Event structure.
@author: Collin
'''

from collections import namedtuple
from operator import *
import itertools as it

class Shape(object):
    """
    Top-level shape that basically only defines some basic
    attributes of shapes
    """
    def __init__(self, start, stop, mid, height):
        if start > 0 and stop > 0 and height > 0:
            self.start = start
            self.stop = stop
            self.mid = mid
            self.height = height
            self.top = mid + height/2.0
            self.bottom = mid - height/2.0

        else:
            raise ValueError("Some values were at or below zero: Start = %s, Stop = %s, Height = %s"%(start, stop, height))


class ExonRectangle(Shape):

    def __init__(self, start, stop, mid, height):
        # call parent init
        super(ExonRectangle, self).__init__(start, stop, mid, height)

        # Declare derived attributes
        Coord = namedtuple('Coord', 'x y')  # makes it easy to access x, y coordinate data
        self.top_left = Coord(x = self.start, y = self.top)
        self.top_right = Coord(x = self.stop, y = self.top)
        self.bottom_left = Coord(x = self.start, y = self.bottom)
        self.bottom_right = Coord(x = self.stop, y = self.bottom)
        self.left = Coord(x = self.start, y = 0)
        self.right = Coord(x = self.stop, y = 0)
        self.exon_start = self.start
        self.exon_stop = self.stop

    def __len__(self):
        """
        allows the len() function to be called
        """
        return int(self.stop - self.start)

    def shift(self, dif):
        """
        The shift method is called when intron down scaling happens and the exon shape needs to move to 'fake'
        coordinates so that the exon looks bigger. The shift method essentially shifts each relevant attribute
        by a parameter passed value.
        """
        self.start += dif
        self.stop += dif
        self.exon_start += dif
        self.exon_stop += dif
        Coord = namedtuple("Coord", "x y")
        self.top_left = Coord(x = self.top_left.x + dif, y = self.top_left.y)
        self.top_right = Coord(x = self.top_right.x + dif, y = self.top_right.y)
        self.bottom_left = Coord(x = self.bottom_left.x + dif, y = self.bottom_left.y)
        self.bottom_right = Coord(x = self.bottom_right.x + dif, y = self.bottom_right.y)
        self.left = Coord(x = self.left.x + dif, y = self.left.y)
        self.right = Coord(x = self.right.x + dif, y = self.right.y)

class JunctionLine(object):
    """
    """
    def __init__(self, exon_list):
        self.exon_lines = []  # straight lines between exons
        self.exon_list = sorted(exon_list, key = lambda x: (x.start, x.stop))

    def createExonLines(self):
        self.exon_lines = []
        for from_exon, to_exon in it.izip(self.exon_list, self.exon_list[1:]):
            self.exon_lines.append([from_exon.right, to_exon.left])

    def get_exon_lines(self):
        return self.exon_lines
