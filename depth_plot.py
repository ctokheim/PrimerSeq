#!/usr/bin/env python
# Copyright (C) 2012-2013  Collin Tokheim
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

# graphing imports
import matplotlib
# matplotlib.use('Agg')  # prevent figure pop-ups
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText   # import anchored text
import itertools as it
import argparse
import os
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, DrawingArea, HPacker  # import for offset text box

import traceback  # debugging import

# important packages
# import numpy as np

# own imports
import wig
import utils


class DropFormatter(ticker.ScalarFormatter):
    """
    git rid of bottom left labels
    """
    def __call__(self, x, pos=None):
        if pos==0: return ''
        return ticker.ScalarFormatter.__call__(self, x, pos=None)


class NoFormatter(ticker.ScalarFormatter):
    """
    remove all tick labels
    """
    def __call__(self, x, pos=None):
        return ''


def offset_text(ax, txt, loc, bbox_to_anchor=(-.03, .65), text_props={"color": "k", "size": "13"}):
    """This function is a copy of the one in draw.py"""
    box1 = TextArea(txt, textprops=text_props)

    box = HPacker(children=[box1],
                  align="center",
                  pad=2, sep=5)

    anchored_box = AnchoredOffsetbox(loc=loc,
                                     child=box, pad=0.,
                                     frameon=False,
                                     prop=dict(size=5),
                                     bbox_to_anchor=bbox_to_anchor,
                                     bbox_transform=ax.transAxes,
                                     borderpad=0.1,
                                     )
    ax.add_artist(anchored_box)


def draw_text(ax, txt):
    """
    Add text in right corner that identifies the data
    """

    at = AnchoredText(txt,
                      loc=1, prop=dict(size=13), frameon=True,)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)


# add commas to an integer for display on figure
def addCommas(myNum):
    """
    A long function for a simple task, add commas to a number
    """
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


def scale_depth(depth_dictionary, start, stop, step):
    list_of_ranges = []
    all_pos = sorted(it.chain(start, stop))
    for i in range(len(all_pos)-1):
        if not i % 2:
            list_of_ranges.append(range(all_pos[i], all_pos[i+1]))
        else:
            list_of_ranges.append(range(all_pos[i], all_pos[i+1], step))

    scaled_dict = {}
    pos_counter = start[0]
    for pos in it.chain.from_iterable(list_of_ranges):
        scaled_dict[pos_counter] = depth_dictionary.get(pos, 0)
        pos_counter += 1

    return scaled_dict, start[0], pos_counter


## Plot generating functions ##
def generate_plot(ax, bw, chr, start, stop, options):
    """
    Calls the matplotlib hist method to create the histogram.
    """
    #depth_counts, max_count = sam.get_depth(bam, chr, start, stop)
    wig_obj = wig.Wig(bw, ext='wig')
    tmp_start, tmp_stop = (start[0], stop[-1]) if type(start) == type(tuple()) else (start, stop)
    wig_obj.extractBigRegion(chr, tmp_start, tmp_stop)
    wig_obj.load_wig_file()
    depth_dict = wig_obj.get_annotation()
    if options['step'] != 1: depth_dict, tmp_start, tmp_stop = scale_depth(depth_dict, start, stop, options['step'])
    max_count = max(depth_dict.itervalues()) if len(list(depth_dict.itervalues())) > 0 else 1
    # depth_counts = list(it.chain.from_iterable([key] * depth_dict[key] for key in depth_dict)) if depth_dict else [0]
    # ax.hist(depth_counts, range(tmp_start, tmp_stop + 1), facecolor='k')
    if max_count > 1:
        ax.hist(map(lambda x: (x[0] + x[1]) / 2., depth_dict.keys()),
                # range(tmp_start, tmp_stop + 1),
                wig_obj.bins,
                weights=depth_dict.values(),
                facecolor='k')
    return max_count, tmp_start, tmp_stop


def read_depth_plot(options):
    if type(options['position']) == type(list()):
        chr = utils.get_chr(options['position'][0])
        start, stop = zip(*map(lambda x: utils.get_pos(x), options['position']))
    else:
        chr = utils.get_chr(options['position'])
        start, stop = utils.get_pos(options['position'])
    bigwigs = options['bigwig'].split(',')
    num_subplots = len(bigwigs)  # num of bam files equals number of subplots
    fig, axes = plt.subplots(num_subplots, 1, sharex=True, sharey=True, figsize=(6, options['size'] * num_subplots))
    gray = (0.9, 0.9, 0.9)

    # iterate over subplots (bigwig files)
    max_count_holder = 0
    if num_subplots == 1:
        # axes.set_title('Read Depth Plot on %s' % chr)
        iterable = [axes]
    else:
        # axes.flat[0].set_title('Read Depth Plot on %s' % chr)
        iterable = axes.flat
    for i, ax in enumerate(iterable):
        #ax.locator_params(nbins=2)
        ax.yaxis.set_label_text('')

        # set bg
        ax.patch.set_facecolor(gray)
        ax.patch.set_edgecolor(gray)
        ax.grid()

        # plot/label
        max_count, real_start, real_stop = generate_plot(ax, bigwigs[i], chr, start, stop, options)  # does the actual work
        draw_text(ax, '%s -- ' % options['gene'] + os.path.splitext(os.path.basename(bigwigs[i]))[0])

        # format options
        ax.xaxis.grid(color='white', linestyle='--', linewidth=1.5)
        ax.yaxis.grid(color='white', linestyle='--', linewidth=1.5)
        ax.xaxis.set_major_formatter(DropFormatter())
        ax.yaxis.set_major_formatter(DropFormatter())
        ax.set_axisbelow(True)

        # hide some ugly lines
        for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
            line.set_color(gray)

        # set y-axis
        if max_count > max_count_holder:
            ax.set_ylim(0, 1.5 * max_count)
            ax.set_yticks([0, int( .375 * max_count ), int( .75 * max_count ), int( 1.125 * max_count ), int(1.5 * max_count)])
            max_count_holder = max_count

        # set x-axis options
        ax.set_xlim(real_start, real_stop)     # set x limits
        ax.set_xticks([real_start, real_stop])   # explicitly set ticks
        ax.xaxis.set_ticklabels(map(addCommas, [real_start, real_stop]))   # make nice looking text for labels
        ax.get_xticklabels()[0].set_horizontalalignment('left')
        ax.get_xticklabels()[1].set_horizontalalignment('right')

        # make text box to display chromosome information
        if i == num_subplots - 1:
            offset_text(ax, '%s:' % chr, 3, (-.15, -.16))

        # adjust spacing between subplots
        fig.subplots_adjust(wspace=0.05, hspace=0.05)

        # save figure
        plt.savefig(options['output'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bigwig', action='store', required=True, help='Comma separated list of bigwig files')
    parser.add_argument('-p', dest='position', action='store', required=True, help='Zero based position chr:start-stop')
    parser.add_argument('-g', dest='gene', action='store', required=True, help='the name')
    parser.add_argument('-s', dest='size', action='store', default=2.)
    parser.add_argument('--step', dest='step', action='store', type=int, default=1)
    parser.add_argument('-o', dest='output', action='store', required=True, help='path to output figure')
    options = vars(parser.parse_args())

    if len(options['position']) > 1:
        options['position'] = options['position'].split(',')

    read_depth_plot(options)
