'''
Created on Feb 8, 2012

@author: Collin
'''
# external dependencies
import numpy as np
import matplotlib
from matplotlib import rc
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import math
from collections import Counter, namedtuple
import traceback
import sys
import argparse
import itertools as it
from operator import *

# custom imports
from shapes import ExonRectangle, JunctionLine


def addCommas(myNum):
    strNum = str(int(myNum))[::-1]
    numLength = len(strNum)

    # figure out the number of commas needed
    numCommas = numLength / 3

    # find the remainder
    remainder = numLength % 3

    # get 3 digit substrings
    subStrings = []
    i = 1
    for i in range(0, numCommas):
        tmp = strNum[i*3:i*3+3]
        revTemp = tmp[::-1]
        subStrings = [revTemp] + subStrings

    if not remainder == 0:
        # add first digits in number
        tmp = strNum[-remainder:]
        revTmp = tmp[::-1]
        subStrings = [revTmp] + subStrings

    # return comma seperated integer
    return ",".join(subStrings)


def scale_intron_length(exonShapes, scale):
    firstFlag = True
    tmpList = []
    previous = None
    previousActual = None
    for index, exonRect in enumerate(exonShapes):
        if exonRect == None:
            tmpList.append(None)
        elif not firstFlag == True:
            dist = exonRect.start - previousActual.stop # distance from the previous (non-modified) exon
            tmpStart = int(dist/scale) + previous.stop # scale distance between exons and plot relative to last exon
            dif = tmpStart - exonRect.start
            previousActual = ExonRectangle(start = exonRect.start, stop = exonRect.stop, mid = 0, height = 1)  # the actual coordinates of the exon
            exonRect.shift(dif)
            tmpList.append(exonRect)
            previous = exonRect
        else:
            previousActual = exonRect
            previous = exonRect
            tmpList.append(exonRect)
            firstFlag = False
    return tmpList  # return exon shape objects accounting for positions changing do to scaling introns


def editAxis(ax, include=[], notInclude=[]):
    for loc, spine in ax.spines.iteritems():
        if loc in include:
            spine.set_position(('outward', 10))  # outward by 10 points
        elif loc in notInclude:
            spine.set_color('none')  # don't draw spine
        else:
            raise ValueError('unknown spine location: %s' % loc)


def exonDrawSubplot(ax, exons):  # exons was coords
    # move/remove axis
    includedAxes = ["bottom"]
    notIncludedAxes = ["left", "right", "top"]
    editAxis(ax, include=includedAxes, notInclude=notIncludedAxes)  # funtion to remove/modify axis

    # turn off ticks where there is no line
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # get ExonRectangle objects
    exonRectangles = []
    for myIndex, (tmpStart, tmpStop) in enumerate(exons):
        tmpRect = ExonRectangle(start=tmpStart, stop=tmpStop, mid=0, height=1)
        exonRectangles.append(tmpRect)  # holds the exons to be drawn
    exonRectangles.sort(key=lambda x: (x.start, x.stop))

    # scale intron length
    # exonRectangles = scale_intron_length(exonRectangles, plotParameters.scale)

    # plot exons
    patches = []
    exonCount = 0
    for tmpExon in exonRectangles:
        art = mpatches.Rectangle(np.array([tmpExon.exon_start, tmpExon.bottom_left.y]), tmpExon.exon_stop - tmpExon.exon_start, tmpExon.height)
        art.set_clip_on(False)
        patches.append(art)
        exonCount += 1

    collection = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1)

    # make color scheme size match the number of exons
    color_scheme = [[0, 0, 0]] * len(exons)

    # plot exons
    collection.set_facecolor(color_scheme)
    ax.add_collection(collection)

    # add lines when applicable
    my_junction_line = JunctionLine(exonRectangles)
    my_junction_line.createExonLines()
    exon_lines = my_junction_line.get_exon_lines()

    # exonLines could be emptry
    for tmp_line in exon_lines:
        ax.plot([tmp_line[0][0], tmp_line[1][0]], [tmp_line[0][1], tmp_line[1][1]], linewidth=1)

    # set axis properties
    firstExonStart = exonRectangles[0].start
    lastExonEnd = exonRectangles[-1].stop

    # find min/max
    minLabel = exons[0].split(":")[1].split("-")[0]
    maxLabel = exons[-1].split(":")[1].split("-")[1]

    ax.set_xlim(firstExonStart, lastExonEnd)
    ax.set_xticks([firstExonStart, lastExonEnd])  # edited
    ax.xaxis.set_ticklabels(["%s" % (addCommas(minLabel)), "%s" % (addCommas(maxLabel))])  # prevents scientific notation and provide scale effect for axis

    ax.set_ylim(-2, 2.25)
    plt.gca().axes.get_yaxis().set_visible(False)

    return ax


def main(tx_paths, counts):
    # configurations
    matplotlib.rcParams['font.size'] = 16  # edit font size of title
    matplotlib.rcParams['xtick.labelsize'] = 13
    matplotlib.rcParams['ytick.labelsize'] = 13
    #matplotlib.rc("lines", linewidth=optionDict["thick"])

    # sort the potential paths by normalized counts
    tmp = zip(counts, tx_paths)
    tmp.sort(key=lambda x: -float(x[0])/(len(x[1]) - 1))
    top_counts, top_tx_paths = zip(*tmp)
    if len(top_tx_paths) >= 5:
        counts = top_counts[:5]
        tx_paths = top_tx_paths[:5]
        num_of_txs = 5
    else:
        counts = top_counts
        tx_paths = top_tx_paths
        num_of_txs = len(top_tx_paths)

    fig, axes = plt.subplots(num_of_txs, 1, sharex=True, sharey=True, figsize=(6, 2 * num_of_txs))

    # loop through and make all subplots in a figure
    for i, ax in enumerate(axes.flat):
        # draw exon structure + junctions
        # exonDrawAxis = exonDrawSubplot(ax, tx_paths[i])
        ### start ###
        exons = tx_paths[i]
    
        # move/remove axis
        includedAxes = ["bottom"]
        notIncludedAxes = ["left", "right", "top"]
        editAxis(ax, include=includedAxes, notInclude=notIncludedAxes)  # funtion to remove/modify axis

        # turn off ticks where there is no line
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        # get ExonRectangle objects
        exonRectangles = []
        for myIndex, (tmpStart, tmpStop) in enumerate(exons):
            tmpRect = ExonRectangle(start=tmpStart, stop=tmpStop, mid=0, height=1)
            exonRectangles.append(tmpRect)  # holds the exons to be drawn
        exonRectangles.sort(key=lambda x: (x.start, x.stop))

        # scale intron length
        # exonRectangles = scale_intron_length(exonRectangles, plotParameters.scale)

        # plot exons
        patches = []
        exonCount = 0
        for tmpExon in exonRectangles:
            art = mpatches.Rectangle(np.array([tmpExon.exon_start, tmpExon.bottom_left.y]), tmpExon.exon_stop - tmpExon.exon_start, tmpExon.height)
            art.set_clip_on(False)
            patches.append(art)
            exonCount += 1

        collection = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1)

        # make color scheme size match the number of exons
        color_scheme = [[0, 0, 0]] * len(exons)

        # plot exons
        collection.set_facecolor(color_scheme)
        ax.add_collection(collection)

        # add lines when applicable
        my_junction_line = JunctionLine(exonRectangles)
        my_junction_line.createExonLines()
        exon_lines = my_junction_line.get_exon_lines()

        # exonLines could be emptry
        for tmp_line in exon_lines:
            ax.plot([tmp_line[0][0], tmp_line[1][0]], [tmp_line[0][1], tmp_line[1][1]], linewidth=1)

        # set axis properties
        firstExonStart = exonRectangles[0].start
        lastExonEnd = exonRectangles[-1].stop

        # find min/max
        minLabel = exons[0][0]
        maxLabel = exons[-1][1]

        ax.set_xlim(firstExonStart, lastExonEnd)
        ax.set_xticks([firstExonStart, lastExonEnd])  # edited
        ax.xaxis.set_ticklabels(["%s" % (addCommas(minLabel)), "%s" % (addCommas(maxLabel))])  # prevents scientific notation and provide scale effect for axis

        ax.set_ylim(-2, 2.25)
        plt.gca().axes.get_yaxis().set_visible(False)
        ### END ###

        editAxis(ax, include=["left"], notInclude=["bottom", "top", "right"])  # removes x-axis spline
        #exonDrawAxis.get_xaxis().set_visible(False)  # set axis invisible
        #exonDrawAxis.set_title(exonCoordinates[exonIndex][0] + " " + gene) # add title

    fig.subplots_adjust(hspace=.05, wspace=.05)  # change subplot spacing
    # plt.show()
    plt.savefig('example.png')


if __name__ == '__main__':
    pass  # at this point no command line call
