# -*- coding: utf-8 -*-


"""
==========================================
Null Models for Robust Functional Networks
==========================================

:Author:
    Moritz Emanuel Beber
:Date:
    2011-03-21
:Copyright:
    Copyright |c| 2011 Jacobs University Bremen gGmbH, all rights reserved.
:File:
    nullmodel.py


.. |c| unicode:: U+A9
"""


import random
import multiprocessing
import numpy
import networkx as nx
import meb.utils.statistics as stats
import meb.utils.network.randomisation as net_rnd
import meb.utils.network.subgraphs as net_subs

from meb.utils.mathfuncs import binomial_coefficient


class NullModel(nx.DiGraph):
    """
    """

    def __init__(self, density, nodes_in=8, nodes_middle=20, nodes_out=8,
            num_modules=0, p_out=1):
        """
        Directed Flow Network Null Model.
        """
        from .classes import ParameterManager
        nx.DiGraph.__init__(self)
        self.parameters = ParameterManager(nodes_in, nodes_middle, nodes_out)
        self.possible_edges =\
                self.parameters.nodes_in * self.parameters.nodes_middle +\
                self.parameters.nodes_middle * (self.parameters.nodes_middle - 1)\
                + self.parameters.nodes_middle * self.parameters.nodes_out
        self.density = float(density)
        self.actual_edges = self.density * float(self.possible_edges)
        self.mtf_counts = None
        self.random_ensemble = None

    def copy(self):
        other = NullModel(self.density)
        other.parameters = self.parameters
        other.possible_edges = self.possible_edges
        other.actual_edges = self.actual_edges
        other.add_nodes_from(n for n in self)
        other.add_edges_from((u, v) for (u, v) in self.edges_iter())
        return other

    def generate_network(self, probability=0.0):
        """
        Use the given parameters to construct a randomly initialised flow
        network.
        """
        self.clear()
        for n in xrange(self.parameters.nodes_total):
            self.add_node(n)
        if probability:
            self._probabilistic_init(probability)
        if self.actual_edges > 0:
            self._fixed_init()

    def _fixed_init(self):
        """
        Initialises the network with a fixed number of edges.
        """
        while self.size() < self.actual_edges:
            (u, v) = random.sample(xrange(self.parameters.nodes_total), 2)
            if not self.has_edge(u, v) and\
                    self.parameters.conserves_structure(u, v):
                self.add_edge(u, v)

    def _probabilistic_init(self, probability):
        """
        Initialises the network with a link probability
        """
        for u in xrange(self.parameters.nodes_end_middle):
            for v in xrange(self.parameters.nodes_in,
                    self.parameters.nodes_total):
                if self.parameters.conserves_structure(u, v) and\
                        random.random() < probability:
                    self.add_edge(u, v)

    def generate_random_ensemble(self, random_number=100, structure=False,
            modularity=False):
        """
        Generates randomised versions of the network that may or may not regard
        certain subsidiary conditions.
        """
        pool = multiprocessing.Pool()
        map = pool.map
        if structure:
            self.random_ensemble = map(_structured_rewiring,
                    [(self, x) for x in xrange(random_number)])
        else:
            self.random_ensemble = map(_normal_rewiring,
                    [(self, x) for x in xrange(random_number)])

    def compute_tsp(self):
        if not self.mtf_counts:
            self.mtf_counts = net_subs.triadic_census(self)
        zscores = numpy.zeros(13)
        for mtf_num in xrange(1, 14):
            tricode = net_subs.num2tricode[mtf_num]
            zscores[mtf_num - 1] = stats.compute_zscore(
                    self.mtf_counts.get(tricode, 0.0),
                    [rnd.mtf_counts.get(tricode, 0.0) for rnd in\
                    self.random_ensemble])
        return zscores


def _normal_rewiring(args):
    network = args[0]
    i = args[1]
    rewire = net_rnd.NetworkRewiring()
    (rnd, success) = rewire.randomise(network)
    rnd.mtf_counts = net_subs.triadic_census(rnd)
    print "random", i, "flip success rate", success
    return rnd

def _structured_rewiring(args):
    network = args[0]
    i = args[1]
    rewire = net_rnd.NetworkRewiring()
    rewire.conditions = check_structured
    (rnd, success) = rewire.randomise(network)
    rnd.mtf_counts = net_subs.triadic_census(rnd)
    print "random", i, "flip success rate", success
    return rnd

def check_structured(graph, first, second):
    """
    Standard rewiring conditions as in original theory.
    """
    # curiously the conditions for switching unidirectional and bidirectional
    # links are the same for just slightly different reasons
    if first == second:
        return False
    # disallow links that defy the network structure
    if not graph.parameters.conserves_structure(first[0], second[1]):
        return False
    if not graph.parameters.conserves_structure(second[0], first[1]):
        return False
    # prevent creation of self-links
    if first[0] == second[1]:
        return False
    if second[0] == first[1]:
        return False
    # check if we would create a parallel edge
    if second[1] in graph[first[0]]:
        return False
    if first[1] in graph[second[0]]:
        return False
    # check if we would create a bidirectional link
    # or cover existing reverse link in double link switching
    if first[0] in graph[second[1]]:
        return False
    if second[0] in graph[first[1]]:
        return False
    return True

def predict_complexity(n_in, n_out, k):
    """
    Compute the predicted mean and variance for a randomly set up output pattern
    with a certain parameter setting.

    Parameters
    ----------
    n_in: int
        Number of input nodes.
    n_out: int
        Number of output nodes.
    k: int
        Number of activated output nodes.
    """
    denominator = binomial_coefficient(n_out, k) * 4.0
    nominator = sum(binomial_coefficient(n_out - k, m) * m for m in range(1, k + 1))
    return numpy.reciprocal(nominator / denominator)

