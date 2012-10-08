import networkx as nx
import numpy as np
import sys

def get_biconnected(G):
    """
    Wrapper arround the networkx biconnected_components function. To find out why the biconnected
    components algorithm is useful for finding constitutive exons check the information section or
    wikipedia.
    """

    G_undirected = G.to_undirected()  # make sure undirected graph for biconnected components
    components = filter(lambda x: len(
        x) > 2, map(list, nx.biconnected_components(G_undirected)))  # filter out trivial dyad biconnected components
    return components


def bellman_ford_longest_path(G, num_nodes, visited, weight='weight'):
    """
    Computes the longest path (most total weight) by only considering unexplained edges. That is weights of any edge already in an isoform
    is set to zero. This function tries to minimize the number of isoforms that could possibly be generated based on novel edges.
    Assumes topologically sorted with source node as first in topological sort. Topological sort version runs in O(n+m) instead of O(nm).

    input:
        G - a networkx DAG
        num_nodes - the number of nodes in the entire graph (not just the subgraph)
        weight - the string that represents edge weight for G
    output:
        path - longest path from first node to last
    """
    if not nx.is_biconnected(G.to_undirected()):
        print "ERROR: component should be biconnected"
        sys.exit(1)

    # initialize variables
    sorted_nodes = sorted(G.nodes())
    d = [float('-inf')] * num_nodes  # was len(G)
    d[sorted_nodes[0]] = 0   # initialize source to have 0 distance
    p = [[]] * num_nodes  # was len(G)
    p[sorted_nodes[0]] = [sorted_nodes[0]]
        # initialize source path to be it's self

    # "edge relax"
    for tail_node in sorted_nodes:
        for head_node in G.successors(tail_node):
            # want longest path of unexplained edges, so all explained edges have zero weight
            edge_weight = G[tail_node][head_node][
                weight] if visited[tail_node][head_node] == 0 else 0

            # larger total weight case
            if d[head_node] < d[tail_node] + edge_weight:
                d[head_node] = d[tail_node] + edge_weight
                p[head_node] = p[tail_node] + [head_node]
            # same total weight case, choose edge with greater weight into head node
            elif d[head_node] == (d[tail_node] + edge_weight) and G[tail_node][head_node][weight] > G[p[head_node][-2]][head_node][weight]:
                d[head_node] = d[tail_node] + edge_weight
                p[head_node] = p[tail_node] + [head_node]

    longest_path = p[sorted_nodes[-1]]
    no_newly_visited = reduce(
        lambda x, y: x and (visited[x][y] == 1), longest_path)

    return p[sorted_nodes[-1]], no_newly_visited


def all_path_lengths(G, component, target):
    # get possible lengths
    inc_length, skip_length = [], []
    sub_graph = nx.subgraph(G, component)
    for path in nx.all_simple_paths(sub_graph.copy(), source=component[0], target=component[-1]):
        if target in path:
            inc_length.append(sum(map(lambda x: x[1] - x[0], path[1:-1])))  # length of everything but target exon and flanking constitutive exons
        else:
            skip_length.append(sum(map(lambda x: x[1] - x[0], path[1:-1])))  # length of everything but target exon and flanking constitutive exons
    return inc_length, skip_length


def read_count_em(bcc_paths, sub_graph):
    # useful convenience dicts
    indexToEdge = {i: e for i, e in enumerate(sub_graph.edges())}
    edgeToIndex = {e: i for i, e in enumerate(sub_graph.edges())}

    # set up count/tx info variables
    num_tx = len(bcc_paths)
    num_edges = len(sub_graph.edges())
    read_counts = np.zeros(num_edges)
    for i, val in enumerate(read_counts):
        u, v = indexToEdge[i]
        read_counts[i] = sub_graph[u][v]['weight']
    total_counts = np.sum(read_counts)

    # set up the uncommited matrix Y
    Y = np.zeros((num_edges, num_tx))
    for tx_index, path in enumerate(bcc_paths):
        for i in range(len(path) - 1):
            try:
                Y[edgeToIndex[(path[i], path[i + 1])]][tx_index] = read_counts[
                    edgeToIndex[(path[i], path[i + 1])]]
            except:
                print '*' * 10, path[i], path[i + 1], tx_index, edgeToIndex
                raise

    # set up p the probability array
    p = np.ones(num_tx) * 1 / num_tx

    # start EM
    epsilon = float('inf')
    THRESHOLD = .0001
    while epsilon > THRESHOLD:
        # E-step
        for i, row in enumerate(Y):
            Y[i] = row * p / row.dot(p) * read_counts[i]

        # M-step
        p_new = np.sum(Y, axis=0) / total_counts

        # convergence variable
        epsilon = np.sum(np.abs(p_new - p))

        p = p_new

    tx_counts = total_counts * p
    return tx_counts


def generate_isoforms(G, tx_paths):
    """
    Take in a digraph and output isoforms of its biconnected components (splice modules)

    input
        G - a directed acyclic graph
        tx_paths - a list of paths from the annotated bed file
    output
        paths - a list of isoforms for splice modules
    """
    paths, module_info, naive = [], [], []
    for component in get_biconnected(G):
        # define subgraph of biconnected nodes
        component_subgraph = nx.subgraph(G, component)
        component_paths = []

        # mark all edges as unvisited
        visited_edges = {}
        for u, v in component_subgraph.edges():
            try:
                visited_edges[u][v] = 0
            except:
                visited_edges[u] = {}
                visited_edges[u][v] = 0
        num_visited_edges = 0  # no edges visited to start with

        # mark visited edges if in known annotation
        for path in tx_paths:
            for i in range(len(path) - 1):
                try:
                    if not visited_edges[path[i]][path[i + 1]]:
                        visited_edges[path[i]][path[i + 1]] = 1
                        num_visited_edges += 1
                except:
                    pass

        # Iterate until atleast one isoform explains each edge
        while num_visited_edges < component_subgraph.number_of_edges():
            # find maximum weight path
            path, no_new_edges = bellman_ford_longest_path(component_subgraph.copy(), num_nodes=G.number_of_nodes(), visited=visited_edges)
            if not no_new_edges:
                component_paths.append(path)

            # bookeeping on visited edges
            for i in range(len(path) - 1):
                if visited_edges[path[i]][path[i + 1]] == 0:
                    visited_edges[path[i]][path[i + 1]] = 1
                    num_visited_edges += 1
                G[path[i]][path[i + 1]]['weight'] = G[path[
                    i]][path[i + 1]]['weight'] / 10.0  # down scale path

        tmp = sorted(component)  # make sure nodes are sorted
        start = tmp[0] if tmp[0] != 0 else tmp[1]  # find min non-dummy node
        stop = tmp[-1] if tmp[-1] != len(G) - 1 else tmp[-2]
                                         # find max non-dummy node
        component_paths += [filter(lambda x: x >= start and x <= stop, p) for p in tx_paths]  # add paths known from transcript annotation
        component_paths = map(list, list(
            set(map(tuple, component_paths))))  # remove redundant paths

        # remove dummy nodes if possible
        if tmp[-1] == len(G) - 1:
            component_subgraph.remove_node(tmp[-1])
        if tmp[0] == 0:
            component_subgraph.remove_node(0)

        count_info = read_count_em(component_paths, component_subgraph)
        original_count = sum([component_subgraph[arc_tail][arc_head]['weight'] for arc_tail, arc_head in component_subgraph.edges()])  # total actual read counts that the subgraph has
        module_info.append([[start, stop], list(count_info), original_count])   # module starts at first node and ends at last node
        paths.append(component_paths)  # add paths

    return paths, module_info
