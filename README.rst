A Python implementation of the interval binning scheme
======================================================

These are some utility functions for working with the interval binning scheme
as used in the `UCSC Genome Browser
<http://genome.cshlp.org/content/12/6/996.full>`_. This scheme can be used to
implement fast overlap-based querying of intervals, essentially mimicking an
`R-tree <https://en.wikipedia.org/wiki/R-tree>`_ index.

Note that some database systems natively support spatial index methods such as
R-trees. See for example the `PostGIS <http://postgis.net/>`_ extension for
PostgreSQL.

Although in principle the method can be used for binning any kind of
intervals, be aware that the largest position supported by this implementation
is 2^29 (which covers the longest human chromosome).


Usage
-----

See the `documentation <https://interval-binning.readthedocs.org/>`_.


Installation
------------

To install the latest release via PyPI using pip::

    pip install interval-binning
