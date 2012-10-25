import algorithms as algs
import logging
import json


class ExonSeek(object):
    '''
    This class handles the finding exons with a psi value threshold.
    '''

    def __init__(self, target, splice_graph):
        self.target = target
        self.graph = splice_graph.get_graph()  # convenience variable (could just use splice_graph)
        if self.target not in self.graph.nodes():
            raise ValueError('The target was not found in the graph (likely a bug in the code)')

        self.strand = splice_graph.strand  # convenience variable
        self.splice_graph = splice_graph
        biconnected_comp = filter(lambda x: target in x, algs.get_biconnected(self.graph))
        self.upstream, self.downstream, self.total_components = None, None, None
        self.psi_upstream, self.psi_target, self.psi_downstream = None, None, None

        self.num_of_biconnected = len(biconnected_comp)
        if len(self.graph.predecessors(self.target)) == 0 or len(self.graph.successors(self.target)) == 0:
            print 'what am i doing here'
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
            raise ValueError('I expect there to be either 0, 1, or 2 biconnected components. Received %s' % self.num_of_biconnected)

    def get_info(self):
        '''
        Return all of the important variables in just one method call.
        '''
        return self.upstream, self.downstream, self.component, self.psi_target, self.psi_upstream, self.psi_downstream

    def find_closest_exon_above_cutoff(self, paths, counts, possible_exons, CUT_OFF=.95):
        """
        Progressively step away from the target exon to find a sufficient constitutive exon
        """
        for exon in possible_exons:
            psi = algs.estimate_psi(exon, paths, counts)
            if psi >= CUT_OFF:
                return exon, psi

    def two_biconnected_case(self):
        print 'two case'
        if self.component[0][-1] == self.target:
            before_component, after_component = self.component
        else:
            after_component, before_component = self.component

        # since there is two components I need two subgraphs/paths. One for
        # before and after the target exon (before/after are defined by
        # chromosome position)
        my_before_subgraph = self.graph.subgraph(before_component)
        before_paths, before_counts = algs.generate_isoforms(my_before_subgraph, self.splice_graph)
        my_before_subgraph = self.graph.subgraph(before_component)
        after_paths, after_counts = algs.generate_isoforms(my_before_subgraph, self.splice_graph)

        if self.strand == '+':
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

    def no_biconnected_case(self):
        # add information to log file
        logging.debug('It appears %s has two imediate flanking constitutive exons' % str(self.target))
        if len(self.graph.successors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and successors')
        if len(self.graph.predecessors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and predecessors')

        # define adjacent exons as flanking constitutive since all three (the
        # target exon, upstream exon, and downstream exon) are constitutive
        self.upstream = self.graph.predecessors(self.target)[0] if self.strand == '+' else self.graph.successors(self.target)[0]
        self.downstream = self.graph.successors(self.target)[0] if self.strand == '+' else self.graph.predecessors(self.target)[0]

        # defining two attributes as the same thing seems silly but in a
        # different case with two biconnected components the two components
        # need to be merged into a single self.total_components
        self.total_components = [self.upstream, self.target, self.downstream]
        self.component = self.total_components

    def one_biconnected_case(self):
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
        index = self.component.index(self.target)
        my_subgraph = self.graph.subgraph(self.component)
        paths, counts = algs.generate_isoforms(my_subgraph, self.splice_graph)
        if self.strand == '+':
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

        with open('tmp/isoforms/test.json', 'w') as handle:
            json.dump({'path': paths, 'counts': list(counts)}, handle, indent=4)  # output path information to tmp file
        self.psi_target = algs.estimate_psi(self.target, paths, counts)

    def first_exon_case(self):
        if len(self.graph.predecessors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and predecessors')
        my_subgraph = self.graph.subgraph(self.component)
        paths, counts = algs.generate_isoforms(my_subgraph, self.splice_graph)
        if self.strand == '+':
            self.upstream = self.graph.predecessors(self.target)[0]
            self.psi_upstream = 1.0  # defined by biconnected component alg as constitutive
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts, self.component[1:])
        else:
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts, self.component[1:])
            self.downstream = self.graph.predecessors(self.target)[0]
            self.psi_downstream = 1.0
        self.psi_target = 1.0

    def last_exon_case(self):
        if len(self.graph.successors(self.target)) > 1:
            logging.debug('Conflict between biconnected components and successors')

        possible_const = self.component[:-1]
        possible_const.reverse()  # reverse the order since closer exons should be looked at first
        my_subgraph = self.graph.subgraph(self.component)
        paths, counts = algs.generate_isoforms(my_subgraph, self.splice_graph)
        if self.strand == '+':
            self.upstream, self.psi_upstream = self.find_closest_exon_above_cutoff(paths,
                                                                                   counts, possible_const)
            self.downstream = self.graph.successors(self.target)[0]
            self.psi_downstream = 1.0
        else:
            self.upstream = self.graph.successors(self.target)[0]
            self.psi_upstream = 1.0
            self.downstream, self.psi_downstream = self.find_closest_exon_above_cutoff(paths,
                                                                                       counts, possible_const)
        self.psi_target = 1.0  # the target is constitutive in this case


