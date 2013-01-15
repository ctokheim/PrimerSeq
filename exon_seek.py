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

'''
File: exon_seek.py
Author: Collin Tokheim
Description: exon_seek.py holds the ExonSeek class which searches
for appropriate flanking "constitutive" exons to place primers on.
'''
import algorithms as algs
import logging
import json
import sequence_interval
import utils


class ExonSeek(object):
    '''
    This class handles the finding exons with a psi value threshold. That is,
    it finds flanking exons to place primers where the exon inclusion level is
    above a user-defined value.
    '''

    def __init__(self, target, splice_graph, ID, cutoff, upstream=None, downstream=None):
        '''
        As the purpose of ExonSeek is to flanking constitutive exons, it is
        necessary to know what needs to be "flanked", the target, and have a
        splice graph representation of gene structure (splice_graph). The ID
        variable is meant to prevent overwriting of files.
        '''
        self.id = ID  # id is to prevent overwriting files in self.save_path_info
        self.cutoff = cutoff
        self.target = target  # (start, end)
        self.upstream, self.downstream = upstream, downstream  # if set, the user specified upstream or downstream exons
        self.graph = splice_graph.get_graph()  # convenience variable (could just use splice_graph)
        if self.target not in self.graph.nodes():
            raise utils.PrimerSeqError('Error: The target was not found in the splice graph')

        self.strand = splice_graph.strand  # convenience variable
        self.splice_graph = splice_graph
        biconnected_comp = filter(lambda x: target in x, algs.get_biconnected(self.graph))
        self.total_components = None  # these will be defined after calling methods
        self.psi_upstream, self.psi_target, self.psi_downstream = None, None, None  # these will be defined after calling methods
        self.all_paths = None

        self.num_of_biconnected = len(biconnected_comp)
        if len(self.graph.predecessors(self.target)) == 0 or len(self.graph.successors(self.target)) == 0:
            self.component = None  # no flanking exon case
        elif self.num_of_biconnected == 0:
            self.no_biconnected_case()
        elif self.num_of_biconnected == 1:
            self.component = sorted(biconnected_comp[0], key=lambda x: (x[0], x[1]))  # make sure component is sorted by position
            self.one_biconnected_case()
        elif self.num_of_biconnected == 2:
            self.component = map(lambda x: sorted(x, key=lambda y: (y[0], y[1])), biconnected_comp)
            self.two_biconnected_case()
        else:
            raise ValueError('Error: There to be either 0, 1, or 2 biconnected components. Received %s' % self.num_of_biconnected)

    def get_info(self):
        '''
        Return all of the important variables in just one method call.
        '''
        return self.all_paths, self.upstream, self.downstream, self.component, self.psi_target, self.psi_upstream, self.psi_downstream

    def save_path_info(self, p, cts):
        '''
        Save information about isoforms and their read counts into a json file.
        '''
        with open('tmp/isoforms/' + str(self.id) + '.json', 'w') as handle:
            json.dump({'path': p, 'counts': list(cts)}, handle, indent=4)  # output path information to tmp file
        self.path, self.counts = p, cts  # store fore programatic access

    def find_closest_exon_above_cutoff(self, paths, counts, possible_exons):
        """
        Progressively step away from the target exon to find a sufficient constitutive exon
        """
        for exon in possible_exons:
            psi = algs.estimate_psi(exon, paths, counts)
            if psi >= self.cutoff:
                return exon, psi

    def two_biconnected_case(self):
        '''
        This is a case where the target exon is constitutive but has two
        flanking biconnected components. Meaning estimating psi for both
        the upstream and downstream exon is necessary
        '''
        print 'two biconnected case'
        if self.component[0][-1] == self.target:
            before_component, after_component = self.component
        else:
            after_component, before_component = self.component

        # since there is two components I need two subgraphs/paths. One for
        # before and after the target exon (before/after are defined by
        # chromosome position)
        before_all_paths = algs.AllPaths(self.splice_graph, before_component, self.target, self.splice_graph.chr)
        before_all_paths.trim_tx_paths()
        before_paths, before_counts = before_all_paths.estimate_counts()
        after_all_paths = algs.AllPaths(self.splice_graph, after_component, self.target, self.splice_graph.chr)
        after_all_paths.trim_tx_paths()
        after_paths, after_counts = after_all_paths.estimate_counts()

        if self.upstream and self.downsteam:
            if self.strand == '+':
                self.psi_upstream = algs.estimate_psi(self.upstream, before_paths, before_counts)
                self.psi_downstream = algs.estimate_psi(self.downstream, after_paths, after_counts)
            elif self.strand == '-':
                self.psi_upstream = algs.estimate_psi(self.upstream, after_paths, after_counts)
                self.psi_downstream = algs.estimate_psi(self.downstream, before_paths, before_counts)
        elif self.strand == '+':
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(before_paths,
                                                                                   before_counts,
                                                                                   list(reversed(before_component[:-1])))
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(after_paths,
                                                                                       after_counts,
                                                                                       after_component[1:])
        else:
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(after_paths,
                                                                                   after_counts,
                                                                                   after_component[1:])
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(before_paths,
                                                                                       before_counts,
                                                                                       list(reversed(before_component[:-1])))
        self.total_components = before_component[:-1] + after_component
        self.psi_target = 1.0

        # handle the combined components
        tmp_start_ix = self.total_components.index(self.upstream) if self.splice_graph.strand == '+' else self.total_components.index(self.downstream)
        tmp_end_ix = self.total_components.index(self.downstream) if self.splice_graph.strand == '+' else self.total_components.index(self.upstream)
        self.all_paths = algs.AllPaths(self.splice_graph, self.total_components[tmp_start_ix:tmp_end_ix], self.target, self.splice_graph.chr)
        self.all_paths.trim_tx_paths()
        self.all_paths.set_all_path_coordinates()
        paths, counts = self.all_paths.estimate_counts()  # used to be self.before_all_paths
        self.save_path_info(paths, counts)

    def no_biconnected_case(self):
        '''
        Case where the target, upstream, and downstream exons are all
        constitutive. Thus just return the immediate upstream and downstream
        exon along with original target.
        '''
        print 'no biconnected case'
        # add information to log file
        logging.debug('It appears %s has two imediate flanking constitutive exons' % str(self.target))
        if len(self.graph.successors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and successors')
        if len(self.graph.predecessors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and predecessors')

        if self.upstream and self.downstream:
            tmp_upstream = self.graph.predecessors(self.target)[0] if self.strand == '+' else self.graph.successors(self.target)[0]

        # define adjacent exons as flanking constitutive since all three (the
        # target exon, upstream exon, and downstream exon) are constitutive
        tmp_upstream = self.graph.predecessors(self.target)[0] if self.strand == '+' else self.graph.successors(self.target)[0]
        tmp_downstream = self.graph.successors(self.target)[0] if self.strand == '+' else self.graph.predecessors(self.target)[0]
        if self.upstream and self.downstream:
            if self.upstream != tmp_upstream or self.downstream != tmp_downstream:
                # raise error if the user defined exon does not match expectation
                raise utils.PrimerSeqError('Error: Flanking exon choice too far from target exon')
        self.upstream = tmp_upstream  # assign upstream after user-defined exon check
        self.downstream = tmp_downstream  # assign downstream after user-defined exon check

        # defining two attributes as the same thing seems silly but in a
        # different case with two biconnected components the two components
        # need to be merged into a single self.total_components
        self.total_components = [self.upstream, self.target, self.downstream]
        self.component = self.total_components

        # create a dummy all paths variable even though there is only one path
        self.all_paths = algs.AllPaths(self.splice_graph, self.component, self.target, self.splice_graph.chr)
        self.all_paths.trim_tx_paths()
        self.all_paths.set_all_path_coordinates()

        # only one isoform, so read counts do not really matter
        paths, counts = self.all_paths.estimate_counts()
        self.save_path_info(paths, counts)

        # since the upstream, target, and downstream exon are constitutive then
        # they all have inclusion of 1.0
        self.psi_target, self.psi_upstream, self.psi_downstream = 1.0, 1.0, 1.0

    def one_biconnected_case(self):
        '''
        Target exon could be cons_titutive or alternatively spliced.
        '''
        if self.target == self.component[0]:
            # constitutive exon of biconnected component, exons with > start pos are
            # not constitutive. However, the immediate preceding exon will be
            # constitutive
            self.first_exon_case()
        elif self.target == self.component[-1]:
            # constitutive exon of biconnected component, exons with < start pos are not
            # constitutive. However, the immediate successor exon will be
            # constitutive.
            self.last_exon_case()
        else:
            # non-constitutive exon case
            self.non_constitutive_case()
        self.total_components = self.component

    def non_constitutive_case(self):
        '''
        In this case, I also estimate the psi for the target exon since
        it is alternatively spliced. Both upstream and downstream exons are
        checked for the closest sufficiently included exon.
        '''
        print 'non-constitutive case'
        index = self.component.index(self.target)

        # get tx path information
        self.all_paths = algs.AllPaths(self.splice_graph, self.component, self.target, self.splice_graph.chr)
        self.all_paths.trim_tx_paths()
        self.all_paths.set_all_path_coordinates()
        paths, counts = self.all_paths.estimate_counts()

        if self.upstream and self.downstream:
            # known flanking exon case
            self.psi_upstream = algs.estimate_psi(self.upstream, paths, counts)
            self.psi_downstream = algs.estimate_psi(self.downstream, paths, counts)
        elif self.strand == '-':
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts,
                                                                                   self.component[index + 1:])
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts,
                                                                                       list(reversed(self.component[:index])))
        else:
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts,
                                                                                   list(reversed(self.component[:index])))
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts,
                                                                                       self.component[index + 1:])
        self.save_path_info(paths, counts)
        self.psi_target = algs.estimate_psi(self.target, paths, counts)

    def first_exon_case(self):
        '''
        Case where the target and one flanking exon is constitutive.
        '''
        print 'first exon case'
        if len(self.graph.predecessors(self.target)) > 1:
            logging.debug('Error: Conflict between biconnected components and predecessors')

        # get tx path information
        self.all_paths = algs.AllPaths(self.splice_graph, self.component, self.target, self.splice_graph.chr)
        self.all_paths.trim_tx_paths()
        self.all_paths.set_all_path_coordinates()
        paths, counts = self.all_paths.estimate_counts()

        if self.upstream and self.downstream:
            # user defined flanking exon case
            if self.strand == '+' and self.graph.predecessors(self.target)[0] == self.upstream:
                self.psi_upstream = 1.0
                self.psi_downsteam = algs.estimate_psi(self.downstream, paths, counts)
            elif self.strand == '-' and self.graph.predecessors(self.target)[0] == self.downstream:
                self.psi_downstream = 1.0
                self.psi_upstream = algs.estimate_psi(self.upstream, paths, counts)
            else:
                raise utils.PrimerSeqError('Error: Flanking exon choice too far from target exon')
        elif self.strand == '+':
            self.upstream = self.graph.predecessors(self.target)[0]
            self.psi_upstream = 1.0  # defined by biconnected component alg as constitutive
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts, self.component[1:])
            self.save_path_info([[self.upstream] + p for p in paths], counts)  # add const. upstream exon to all paths
        else:
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts, self.component[1:])
            self.downstream = self.graph.predecessors(self.target)[0]
            self.psi_downstream = 1.0
            self.save_path_info([p + [self.downstream] for p in paths], counts)  # add const. downstream exon to all paths
        self.psi_target = 1.0

    def last_exon_case(self):
        '''
        Case where the target and one flanking exon are constitutive.
        '''
        print 'last exon case'
        if len(self.graph.successors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and successors')

        possible_const = self.component[:-1]
        possible_const.reverse()  # reverse the order since closer exons should be looked at first

        # get tx path information
        self.all_paths = algs.AllPaths(self.splice_graph, self.component, self.target, self.splice_graph.chr)
        self.all_paths.trim_tx_paths()
        self.all_paths.set_all_path_coordinates()
        paths, counts = self.all_paths.estimate_counts()

        if self.upstream and self.downstream:
            # user defined flanking exon case
            if self.strand == '+':
                self.psi_upstream = algs.estimate_psi(self.upstream, paths, counts)
                self.psi_downstream = 1.0
            elif self.strand == '-':
                self.psi_upstream = 1.0
                self.psi_downstream = algs.estimate_psi(self.downstream, paths, counts)
        if self.strand == '+':
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts, possible_const)
            self.downstream = self.graph.successors(self.target)[0]
            self.psi_downstream = 1.0
            self.save_path_info([p + [self.downstream] for p in paths], counts)  # add const. downstream exon to all paths
        else:
            self.upstream = self.graph.successors(self.target)[0]
            self.psi_upstream = 1.0
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts, possible_const)
            self.save_path_info([[self.upstream] + p for p in paths], counts)  # add const. upstream exon to all paths
        self.psi_target = 1.0  # the target is constitutive in this case
