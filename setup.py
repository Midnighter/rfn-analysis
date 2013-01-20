#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
RFN Analysis Package
====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-01-20
:Copyright:
    Copyright |c| 2013 Jacobs University Bremen gGmbH, all rights reserved.
:File:
    setup.py
"""


from distutils.core import setup


def get_version():
    import sys
    sys.path.insert(0, __file__)
    ra = __import__(name="rfn_analysis")
    return ra.__version__


setup(
    name = "rfn-analysis",
    version = get_version(),
    description = "flow network analysis package",
    author = "Moritz Emanuel Beber",
    author_email = "m (dot) beber (at) jacobs (dash) university (dot) de",
    url = "https://github.com/Midnighter/rfn-analysis",
    packages = [
            "rfn_analysis"
            ],
    )

