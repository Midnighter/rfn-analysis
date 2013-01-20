# -*- coding: UTF-8 -*-


"""
=================================
File Assembly and Data Extraction
=================================

:Author:
    Moritz Emanuel Beber
:Date:
    2010-09-13
:Copyright:
    Copyright |c| 2010 Jacobs University Bremen gGmbH, all rights reserved.
:File:
    extraction.py


.. |c| unicode:: U+A9
"""


__all__ = ["extract_all", "find_files"]


import os
import random
import numpy
import networkx as nx

from operator import attrgetter
# import necessary for unpickling
#from . import classes as rfn_cls


def find_files(directory, file_pattern, sample=0):
    """
    Given a directory collect filenames matching a pattern.

    Parameters
    ----------
    directory: str
        A directory path.
    file_pattern: regex
        A compiled regular expression from the re module.
    sample: int (optional)
        Return only a random sample of networks.
    """
    files = os.listdir(directory)
    net_files = [os.path.join(directory, filename) for filename in files
            if file_pattern.match(filename)]
    if sample > 0:
        net_files = random.sample(net_files, sample)
    return net_files

def extract_all((networks, net_type, setup, args)):
    """
    Open a pickled network and extract data from it.

    Parameters
    ----------
    networks: iterable
        Iterable of filenames
    net_type: str
        Specific string labelling this data
    setup: str
        String describing the parameters
    args: tuple
        Attributes to be collected (e.g., "robustness")
    """
    getter = attrgetter(*args)
    z_getter = attrgetter("zscores")
    res = list()
    for filename in networks:
        try:
            net = nx.read_gpickle(filename)
        except (IOError, EOFError):
            print "failed to load network file '%s'" % filename
            os.rename(filename, filename + ".failed")
            continue
        results = list(z_getter(net))
        results.extend(list(getter(net)))
        results.append(numpy.mean(net.shortest_paths))
        # stripping .pkl file extension
        results.append(os.path.basename(filename)[:-4])
        results.append(net_type)
        results.append(setup)
        res.append(results)
    return res

