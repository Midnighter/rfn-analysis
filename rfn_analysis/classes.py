# -*- coding: utf-8 -*-


"""
======================================
Classes for Robust Functional Networks
======================================

:Author:
    Moritz Emanuel Beber
:Date:
    2010-09-13
:Copyright:
    Copyright |c| 2010 Jacobs University Bremen gGmbH, all rights reserved.
:File:
    classes.py


.. |c| unicode:: U+A9
"""


__all__ = ["ParameterManager", "RobustFunctionalNetwork"]


import os
import struct
import numpy
import numpy.ma
import warnings
import networkx as nx
import scipy.linalg as la
import meb.utils.statistics as stats
import meb.utils.network.randomisation as net_rnd
import meb.utils.network.subgraphs as net_subs
import meb.utils.network.algorithms as net_algs
import meb.utils.network.community as net_com
from . import nullmodel as null


def print_nets(arg, dirname, fnames):
    """
    use with os.path.walk
    """
    for filename in fnames:
        if filename.endswith("bz2"):
            fname = os.path.join(dirname, filename)
            g = nx.read_gpickle(fname)
            fname = fname.split(".", 1)[0]
            prev = g.draw(fname + "_spectral.pdf", spectral_partition=True)
            g.draw(fname + "_louvain.pdf", prev, louvain_partition=True)


class ParameterManager(object):
    """
    Class that encompasses the general parameters of the robust
    functional networks.
    """

    def __init__(self, nodes_in=8, nodes_middle=20, nodes_out=8, activated_k=4,
            filename=None, *args, **kwargs):
        """
        """
        object.__init__(self, *args, **kwargs)
        self.nodes_in = int(nodes_in)
        self.nodes_middle = int(nodes_middle)
        self.nodes_out = int(nodes_out)
        self.activated_k = int(activated_k)
        if filename:
            with open(filename, "rb") as file_h:
                (self.nodes_in, self.nodes_middle, self.nodes_out,\
                        self.activated_k) = struct.unpack_from("<4H",\
                        file_h.read(), struct.calcsize("<H"))
        self.nodes_end_middle = self.nodes_in + self.nodes_middle
        self.nodes_total = self.nodes_end_middle + self.nodes_out

    def is_input(self, node):
        node = int(node)
        return (node < self.nodes_in)

    def is_middle(self, node):
        node = int(node)
        return (node >= self.nodes_in and node < self.nodes_end_middle)

    def is_output(self, node):
        node = int(node)
        return (node >= self.nodes_end_middle and node < self.nodes_total)

    def conserves_structure(self, u, v):
        """
        Checks that a link u -> v preserves the flow network structure.
        """
        if u == v:
            return False
        elif self.is_input(u) and self.is_input(v):
            return False
        elif self.is_input(u) and self.is_output(v):
            return False
        elif self.is_middle(u) and self.is_input(v):
            return False
        elif self.is_output(u) and self.is_input(v):
            return False
        elif self.is_output(u) and self.is_middle(v):
            return False
        else:
            return True


class RobustFunctionalNetwork(nx.DiGraph):
    """
    RFN network in python using networkx.

    The method load is meant to read the binary format coming from the C++
    simulations. These networks themselves are meant to be saved and loaded
    using the networkx methods, i.e., nx.write_gpickle and nx.read_gpickle.
    """

    def __init__(self, name, parameters, *args, **kw_args):
        """

        """
        nx.DiGraph.__init__(self, name=os.path.abspath(name), *args, **kw_args)
        self.parameters = parameters
        self.flow_error = 0.0
        self.robustness = 0.0
        self.ideal_pattern = None
        self.scalar_complexity = None
        self.binary_complexity = None
        self.iteration = None
        self.possible_edges = None
        self.density = None
        self.initial_connectivity = None
        self.essentiality = None
        self.spectral_modularity = None
        self.spectral_partition = None
        self.louvain_modularity = None
        self.louvain_partition = None
        self.degree_correlation = None
        self.pattern_variance = None
        self.pattern_rank = None
        self.binary_rank = None
        self.members = None
        self.middle_members = None
        self.mtf_counts = None
        self.zscores = None
        self.random_graphs = None
        self.middle_modules = None
        self.mean_overlap = None
        self.shortest_paths = None

    def __str__(self):
        return os.path.basename(self.name)

    def copy(self):
        other = RobustFunctionalNetwork(self.name + " copy", self.parameters)
        other.add_nodes_from(n for n in self)
        other.add_edges_from((u, v) for (u, v) in self.edges_iter())
        return other

    def load(self):
        """
        """
        with open(self.name, "rb") as file_h:
            content = file_h.read()
        # check for failed files
        if len(content) == 0:
            os.remove(self.name)
            raise IOError("empty file")
        # shorten buffer, I think that's more efficient than using offset to
        # search in buffer
        (self.iteration,) = struct.unpack_from("<Q", content)
        content = content[struct.calcsize("<Q"):]
        # skipping seed
        content = content[struct.calcsize("<I"):]
        (self.initial_connectivity,) = struct.unpack_from("<d", content)
        content = content[struct.calcsize("<d"):]
        (self.flow_error, num_edges) = struct.unpack_from("<dI", content)
        content = content[struct.calcsize("<dI"):]
        # now add all edges to network
        edge = struct.Struct("<2I")
        for i in xrange(num_edges):
            (src, tar) = edge.unpack_from(content)
            content = content[edge.size:]
            self.add_edge(src, tar)
        # load ideal pattern
        self.ideal_pattern = numpy.zeros(shape=(self.parameters.nodes_out,
            self.parameters.nodes_in))
        fraction = struct.Struct("<d")
        for i in xrange(self.parameters.nodes_out):
            for j in xrange(self.parameters.nodes_in):
                (value,) = fraction.unpack_from(content)
                self.ideal_pattern[i, j] = value
                content = content[fraction.size:]
        # load robustness if available
        if len(content) == 8:
            (self.robustness,) = struct.unpack_from("<d", content)
            content = content[struct.calcsize("<d"):]
        assert len(content) == 0, "Reading buffer messed up!"

    def analysis(self):
        self.compute_density()
        self.degree_correlation = nx.degree_assortativity_coefficient(self)
        self.compute_complexity()
        self.compute_paths()
        self.compute_overlap()
        self.compute_variance()
        self.pattern_rank = numpy.linalg.matrix_rank(self.ideal_pattern)
        self.binary_rank = numpy.linalg.matrix_rank(self.ideal_pattern > 0)
        try:
            self.compute_modularity()
        except nx.NetworkXError:
            pass
        self.generate_random_ensemble()
        self.compute_zscores()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.compute_essentialities()

    def compute_density(self):
        self.possible_edges =\
                self.parameters.nodes_in * self.parameters.nodes_middle +\
                self.parameters.nodes_middle * (self.parameters.nodes_middle - 1)\
                + self.parameters.nodes_middle * self.parameters.nodes_out
        self.density = float(self.size()) / float(self.possible_edges)

    def compute_modularity(self):
        # use only directed modularity methods
        (self.spectral_modularity, self.spectral_partition) =\
                net_algs.directed_spectral_community_detection(self)
        undir = self.to_undirected()
        self.louvain_partition = net_com.best_partition(undir)
        self.louvain_modularity = net_com.modularity(self.louvain_partition, undir)

    def compute_overlap(self):
        self.middle_modules = None
        self.mean_overlap = None
        if self.size() == 0:
            return
        self.middle_modules = list()
        for j in range(self.parameters.nodes_in):
            if not j in self:
                continue
            middle = set()
            for (u, v) in nx.bfs_edges(self, j):
                if self.parameters.is_middle(u):
                    middle.add(u)
                if self.parameters.is_middle(v):
                    middle.add(v)
            self.middle_modules.append(middle)
        num_in = len(self.middle_modules)
        if num_in == 0:
            return
        # compute overlap
        combinations = (num_in * (num_in - 1)) / 2
        middle_overlap = numpy.zeros(combinations)
        i = 0
        for j in range(num_in - 1):
            for k in range(j + 1, num_in):
                nominator = len(self.middle_modules[j].intersection(self.middle_modules[k]))
                denominator = float(len(self.middle_modules[j].union(self.middle_modules[k])))
                if denominator == 0.0:
                    middle_overlap[i] = numpy.nan
                else:
                    middle_overlap[i] = nominator / denominator
                i += 1
        middle_overlap = numpy.ma.masked_invalid(middle_overlap)
        self.mean_overlap = numpy.mean(middle_overlap[~middle_overlap.mask])

    def compute_paths(self):
        self.shortest_paths = list()
        for src in xrange(self.parameters.nodes_in):
            try:
                (pred, dist) = nx.bellman_ford(self, src)
            except:
                continue
            for tar in xrange(self.parameters.nodes_end_middle,\
                    self.parameters.nodes_total):
                if dist.has_key(tar):
                    self.shortest_paths.append(dist[tar])

    def _get_sle(self):
        """

        """
        sle = numpy.eye(self.parameters.nodes_total,\
            self.parameters.nodes_total)
        for (src, tar) in self.edges_iter():
            sle[tar, src] = -1.0 / float(self.out_degree(src))
        return sle

    def _compute_essentiality(self, original, changed):
        """

        """
        changed = original - changed
        essen = abs(changed).sum()
        essen /= float(2 * self.parameters.nodes_in)
        return essen

    def compute_essentialities(self):
        """

        """
        self.essentiality = dict()
        sle = self._get_sle()
        # get the LU factorised system and the pivoting as a tuple
        original = la.lu_factor(a=sle, overwrite_a=True)
        # input & solution vector
        x = numpy.empty(self.parameters.nodes_total)
        # output pattern
        output_pattern = numpy.zeros(shape=(self.parameters.nodes_out,\
            self.parameters.nodes_in))
        for i in xrange(self.parameters.nodes_in):
            x.fill(0.0)
            x[i] = 1.0
            x = la.lu_solve(original, x, overwrite_b=True)
            output_pattern[:, i] = x[self.parameters.nodes_end_middle:]
        # now compute the change in the output pattern upon deletion of a node
        # or link
        changed_pattern = numpy.empty(shape=(self.parameters.nodes_out,\
            self.parameters.nodes_in))
        # bare copy of the network
        changed_network = self.copy()
        # node essentialities
        for n in self.nodes_iter():
            # store neighbourhood
            links = [(n, succ) for succ in self.adj[n].iterkeys()]
            links.extend([(pred, n) for pred in self.pred[n].iterkeys()])
            # clear the current node
            changed_network.remove_node(n)
            sle = changed_network._get_sle()
            # LU factorise the new system of linear equations
            changed = la.lu_factor(a=sle, overwrite_a=True)
            changed_pattern.fill(0.0)
            for i in xrange(self.parameters.nodes_in):
                x.fill(0.0)
                x[i] = 1.0
                x = la.lu_solve(changed, x, overwrite_b=True)
                changed_pattern[:, i] = x[self.parameters.nodes_end_middle:]
            # compute difference
            self.essentiality[n] = self._compute_essentiality(output_pattern,\
                changed_pattern)
            changed_network.add_edges_from(links)
        # link essentialities
        for (src, tar) in self.edges_iter():
            # clear the current link
            changed_network.remove_edge(src, tar)
            sle = changed_network._get_sle()
            # LU factorise the new system of linear equations
            changed = la.lu_factor(a=sle, overwrite_a=True)
            changed_pattern.fill(0.0)
            for i in xrange(self.parameters.nodes_in):
                x.fill(0.0)
                x[i] = 1.0
                x = la.lu_solve(changed, x, overwrite_b=True)
                changed_pattern[:, i] = x[self.parameters.nodes_end_middle:]
            # compute difference
            self.essentiality[(src, tar)] = self._compute_essentiality(\
                output_pattern, changed_pattern)
            changed_network.add_edge(src, tar)

    def generate_random_ensemble(self, rnd_num=100):
        """
        """
        self.random_graphs = list()
        rewire = net_rnd.NetworkRewiring()
        rewire.conditions = null.check_structured
        self.mtf_counts = net_subs.triadic_census(self)
        for i in xrange(rnd_num):
            (rnd, success) = rewire.randomise(self)
            rnd.mtf_counts = net_subs.triadic_census(rnd)
#            print "random", i, "flip success rate", success
            self.random_graphs.append(rnd)

    def compute_zscores(self):
        self.zscores = numpy.zeros(13)
        for mtf_num in xrange(1, 14):
            tricode = net_subs.num2tricode[mtf_num]
            self.zscores[mtf_num - 1] = stats.compute_zscore(
                    self.mtf_counts.get(tricode, 0.0),
                    [rnd.mtf_counts.get(tricode, 0.0) for rnd in\
                    self.random_graphs])

    def compute_complexity(self):
        """
        Complexity of the output pattern is defined as the sum of the inner product of
        all combinations of output vectors.
        """
        self.scalar_complexity = 0.0
        self.binary_complexity = 0.0
        binary_matrix = (self.ideal_pattern > 0).astype(int)
        for i in xrange(self.parameters.nodes_in - 1):
            for j in xrange(i + 1, self.parameters.nodes_in):
                self.scalar_complexity += numpy.inner(\
                        self.ideal_pattern[:, i],\
                        self.ideal_pattern[:, j])
                self.binary_complexity += numpy.inner(\
                        binary_matrix[:, i],\
                        binary_matrix[:, j])
        # normalise by the number of vectors products, i.e.,
        # nodes_in * (nodes_in - 1) / 2
        norm = float(self.parameters.nodes_in *\
                (self.parameters.nodes_in - 1)) / 2.0
        self.scalar_complexity /= norm
        self.binary_complexity /= (norm * self.parameters.activated_k)
        self.scalar_complexity = numpy.reciprocal(self.scalar_complexity)
        self.binary_complexity = numpy.reciprocal(self.binary_complexity)

    def compute_variance(self):
        """
        Compute the variance of the non-zero output pattern from 1/k.
        """
        self.pattern_variance = numpy.sqrt(numpy.power(
            numpy.reciprocal(self.parameters.activated_k) -\
            self.ideal_pattern[self.ideal_pattern > 0], 2.0).sum())

    def draw_modules(self, filename, previous=None, spectral_partition=False, louvain_partition=False):
        """
        Draw the network using pygraphviz
        """
        try:
            import pygraphviz as pgv
        except ImportError:
            return
        # colours for partitions
        if spectral_partition or louvain_partition:
            # rainbow palette
#            colour_vec = ["#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
#                    "#00FFFFFF", "#0040FFFF", "#8000FFFF"]
            # brewer palette
            colour_vec = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                    "#FFFF33", "#A65628", "#F781BF", "#999999"]
        if previous:
            # only update
            net = pgv.AGraph(previous, name=filename)
            net.has_layout = True
            for node in self.nodes_iter():
                if spectral_partition:
                    for (i, community) in enumerate(self.spectral_partition):
                        if node in community:
                            colour = colour_vec[i]
                elif louvain_partition:
                    colour = colour_vec[self.louvain_partition[node]]
                if spectral_partition or louvain_partition:
                    n = net.get_node(node)
                    n.attr["fillcolor"] = colour
            net.draw(filename)
            return net.to_string()
        else:
            net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr = dict()
        node_attr["style"] = "filled"
        link_attr = dict()
        # add compound nodes
        for node in self.nodes_iter():
            if spectral_partition:
                for (i, community) in enumerate(self.spectral_partition):
                    if node in community:
                        colour = colour_vec[i]
            elif louvain_partition:
                colour = colour_vec[self.louvain_partition[node]]
            if self.parameters.is_input(node):
                label = "I %d" % node
                shape = "invtriangle"
            elif self.parameters.is_middle(node):
                label = "M %d" % node
                shape = "circle"
            else:
                label = "O %d" % node
                shape = "triangle"
            if spectral_partition or louvain_partition:
                net.add_node(node, label=label, fillcolor=colour, shape=shape, **node_attr)
            else:
                net.add_node(node, label=label, shape=shape, **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(u, v, **link_attr)
        net.layout(prog="dot", args="")
        net.draw(filename)
        return net.to_string()

    def draw_layers(self, filename, previous=None, prog="dot", args=""):
        """
        Draw the network using pygraphviz
        """
        try:
            import pygraphviz as pgv
        except ImportError:
            return
        # colours for partitions
        # brewer palette
        colour_vec = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999"]
        if previous:
            # only update
            pass
            net = pgv.AGraph(previous, name=filename)
#            net.has_layout = True
        else:
            net = pgv.AGraph(directed=True, name=filename, strict=True,
                    rankdir="TB")
        node_attr = dict()
        node_attr["style"] = "filled"
        link_attr = dict()
        input_layer = list()
        middle_layer = list()
        output_layer = list()
        # add compound nodes
        for node in self.nodes_iter():
            if self.parameters.is_input(node):
                label = "I %d" % node
                shape = "invtriangle"
                colour = colour_vec[0]
                input_layer.append(node)
            elif self.parameters.is_middle(node):
                label = "M %d" % node
                shape = "circle"
                colour = colour_vec[1]
                middle_layer.append(node)
            else:
                label = "O %d" % node
                shape = "triangle"
                colour = colour_vec[2]
                output_layer.append(node)
#            net.add_node(node, label=label, fillcolor=colour, shape=shape, **node_attr)
            net.add_node(node, label="", fillcolor=colour, shape=shape, **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(u, v, **link_attr)
        # layers defined by subgraphs
        sub_attr = dict()
        net.add_subgraph(input_layer, name="input", rank="source", **sub_attr)
        net.add_subgraph(middle_layer, name="middle", **sub_attr)
        net.add_subgraph(output_layer, name="output", rank="sink", **sub_attr)
        net.layout(prog=prog, args=args)
        net.draw(filename)
        return net.to_string()

