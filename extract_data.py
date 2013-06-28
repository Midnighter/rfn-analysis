#! /usr/bin/env python
# -*- coding: UTF-8 -*-


"""
================================================
Gather Simulation Data for analysis and plotting
================================================

:Author:
    Moritz Emanuel Beber
:Date:
    2011-06-13
:Copyright:
    Copyright |c| 2011 Jacobs University Bremen gGmbH, all rights reserved.
:File:
    extract_data.py


.. |c| unicode:: U+A9
"""


import os
import re
import multiprocessing
import tables

import rfn_analysis as ra


class NetworkStats(tables.IsDescription):
    """
    Row descriptor for pytables.
    """
    mtf_1                = tables.Float64Col()
    mtf_2                = tables.Float64Col()
    mtf_3                = tables.Float64Col()
    mtf_4                = tables.Float64Col()
    mtf_5                = tables.Float64Col()
    mtf_6                = tables.Float64Col()
    mtf_7                = tables.Float64Col()
    mtf_8                = tables.Float64Col()
    mtf_9                = tables.Float64Col()
    mtf_10               = tables.Float64Col()
    mtf_11               = tables.Float64Col()
    mtf_12               = tables.Float64Col()
    mtf_13               = tables.Float64Col()
    flow_error           = tables.Float64Col()
    robustness           = tables.Float64Col()
    scalar_complexity    = tables.Float64Col()
    binary_complexity    = tables.Float64Col()
    iteration            = tables.UInt32Col()
    density              = tables.Float64Col()
    initial_connectivity = tables.Float64Col()
    spectral_modularity  = tables.Float64Col()
    louvain_modularity   = tables.Float64Col()
    degree_correlation   = tables.Float64Col()
    mean_overlap         = tables.Float64Col()
    pattern_variance     = tables.Float64Col()
    pattern_rank         = tables.Float64Col()
    shortest_paths       = tables.Float64Col()
    name                 = tables.StringCol(22)
    type                 = tables.StringCol(12)
    setup                = tables.StringCol(20)


def extract_1d(base, file_pattern, net_stats, setup, header, *attr):
    """
    Multiprocessed network attribute extraction.

    Parameters
    ----------
    base: str
        Base directory path for network files.
    file_pattern: regex
        Compiled regular expression from the re module that filenames must
        match.
    net_stats: tables.Table.row
        A table row that data can be appended to.
    setup: str
        A description of starting parameters.
    *attr:
        Any attributes that should be extracted from the networks.
    """
    pool = multiprocessing.Pool(3)
    arguments = list()
    files = ra.find_files(os.path.join(base, "link_robust"), file_pattern)
    arguments.append((files, "Link Robust", setup, attr))
    files = ra.find_files(os.path.join(base, "node_robust"), file_pattern)
    arguments.append((files, "Node Robust", setup, attr))
    files = ra.find_files(os.path.join(base, "noise_robust"), file_pattern)
    arguments.append((files, "Noise Robust", setup, attr))
    results = pool.map(ra.extract_all, arguments)
    for res in results:
        for stats in res:
            for (i, name) in enumerate(header):
                net_stats[name] = stats[i]
            net_stats.append()

def all_simple_data(source, dest, time, setup):
    """
    Manage the HDF5 file and extraction of data from various sub-directories.

    Parameters
    ----------
    source: str
        A directory path.
    dest: str
        A directory path.
    time: str
        One of: 'flow', 'final', or 'timeline'.
    setup: str
        One of: 'setup' or 'complexity'.
    """

    if time == "flow":
        # networks at threshold
        file_pattern = re.compile(r"sim\d+_flow\.pkl")
    elif time == "final":
        # networks at end of evolution
        file_pattern = re.compile(r"sim\d+_final\.pkl")
    elif time == "timeline":
        # timeline of network evolution
        file_pattern = re.compile(r"sim\d+_\d+\.pkl")
    else:
        raise StandardError("unknown time parameter")

    # attributes from networks to extract additionally
    attr = ("flow_error", "robustness", "scalar_complexity",
            "binary_complexity", "iteration", "density","initial_connectivity",
            "spectral_modularity", "louvain_modularity", "degree_correlation",
            "mean_overlap", "pattern_variance", "pattern_rank")

    header = ["mtf_%d" % i for i in range(1, 14)]
    header.extend(list(attr))
    header.extend(["shortest_paths","name","type","setup"])


    if not os.path.exists(dest):
        os.makedirs(dest)
    h5_file = tables.openFile(os.path.join(dest, "attributes.h5"), mode="a",
            title="Flow Networks Attribute Database")

    if setup == "setups":
        try:
            group = h5_file.root.setups
        except tables.NoSuchNodeError:
            group = h5_file.createGroup("/", "setups", "Variations in the number"\
                    " of input nodes and activated output nodes.")
        if time == "flow":
            try:
                table = group.phase_1
            except tables.NoSuchNodeError:
                table = h5_file.createTable(group, "phase_1", NetworkStats,
                        "Results of evolution phase 1, i.e., pattern recognition.")
        elif time == "final":
            try:
                table = group.phase_2
            except tables.NoSuchNodeError:
                table = h5_file.createTable(group, "phase_2", NetworkStats,
                        "Results of evolution phase 2, i.e., robustness.")

        sub = "standard"
        setup = "Standard"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "4_input"
        setup = "4 Input"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "6_input"
        setup = "6 Input"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "12_input"
        setup = "12 Input"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "10_input"
        setup = "10 Input"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "2_activated"
        setup = "2 Activated"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "6_activated"
        setup = "6 Activated"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "8_activated"
        setup = "8 Activated"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()
    elif setup == "complexity":
        try:
            group = h5_file.root.complexity
        except tables.NoSuchNodeError:
            group = h5_file.createGroup("/", "complexity",
                    "Artificial Ideal Output Patterns")
        if time == "flow":
            try:
                table = group.phase_1
            except tables.NoSuchNodeError:
                table = h5_file.createTable(group, "phase_1", NetworkStats,
                        "Results of evolution phase 1, i.e., pattern recognition.")
        elif time == "final":
            try:
                table = group.phase_2
            except tables.NoSuchNodeError:
                table = h5_file.createTable(group, "phase_2", NetworkStats,
                        "Results of evolution phase 2, i.e., robustness.")
        sub = "equal_complexity"
        setup = "Equal Complexity"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "equal_spread"
        setup = "Equal Spread"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "high_complexity"
        setup = "High Complexity"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()

        sub = "low_complexity"
        setup = "Low Complexity"
        extract_1d(os.path.join(source, sub), file_pattern, table.row, setup,
                header, *attr)
        table.flush()
    else:
        raise StandardError("unknown setup parameter")

    h5_file.close()


if __name__ == "__main__":
    # define directories for in- and output
    base = ""
    source = os.path.join(base, "in")
    dest = os.path.join(base, "out")
    all_simple_data(source, dest, "flow", "complexity")

