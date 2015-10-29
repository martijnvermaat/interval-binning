A Python implementation of the interval binning scheme
======================================================

These are some utility functions for working with the interval binning scheme
as used in the `UCSC Genome Browser
<http://genome.cshlp.org/content/12/6/996.full>`_. This scheme can be used to
implement fast overlap-based querying of intervals, essentially mimicking an
`R-tree <https://en.wikipedia.org/wiki/R-tree>`_ index.

.. note:: Some database systems natively support spatial index methods such as
          R-trees. See for example the `PostGIS <http://postgis.net/>`_
          extension for PostgreSQL.

Although in principle the method can be used for binning any kind of
intervals, be aware that the largest position supported by this implementation
is :math:`2^{29}` (which covers the longest human chromosome).


Usage
-----

Let's say you have a set of intervals `I` in a database system without support
for spatial indexing. Querying `I` on overlap with an interval `q` can be done
as:

.. math::

   \{ i \in I ~|~ \mathrm{overlapping}(i, q) \}

   \text{where}

   \mathrm{overlapping}(i, q) = i.start < q.stop \wedge i.stop > q.start

But this will be slow, even with normal B-tree indexes on `start` and `stop`.

If for each interval `i`, we also store its bin as given by
:func:`.assign_bin` (and we index it), we can get the same result much
faster by pre-filtering on :func:`.overlapping_bins`:

.. math::

   \{ i \in I ~|~ i.bin \in \mathrm{overlapping\_bins}(q) \wedge \mathrm{overlapping}(i, q) \}

Similarly, if `i` must completely contain `q` (or vice versa), you can use
:func:`.containing_bins` (or :func:`.contained_bins`).


SQLAlchemy example
^^^^^^^^^^^^^^^^^^

As a more concrete example, let's consider the following `SQLAlchemy
<http://www.sqlalchemy.org/>`_ model definition for storing gene locations::

    class Gene(Base):
        name = Column(String, primary_key=True)
        start = Column(Integer)
        stop = Column(Integer)
        bin = column(Integer, index=True)

        def __init__(self, name, start, stop):
            self.name = name
            self.start = start
            self.stop = stop
            self.bin = binning.assign_bin(start, stop)

The `bin` column is populated using :func:`.assign_bin` and has an
index. We can ask for all genes spanning the 50,000-50,500 interval::

    >>> session.query(Gene).filter(Gene.start < 50000,
    ...                            Gene.stop > 50500).count()
    78

But that query will be slow, yielding a sequential table scan. Adding indexes
on `Gene.start` or `Gene.stop` will not help much.

The :func:`.containing_bins` function gives us all the bins potentially
containing genes spanning some interval. We can use that to filter on
`Gene.bin` first::

    >>> bins = binning.containing_bins(50000, 50500)
    >>> session.query(Gene).filter(Gene.bin.in_(bins),
    ...                            Gene.start < 50000,
    ...                            Gene.stop > 50500).count()
    78

This query will be much faster because it can use the index on `Gene.bin`. The
filter on bin only gives us a crude pre-selection though, so we still have to
apply the `Gene.start` and `Gene.stop` filters on the (relatively small)
resulting set of genes to get the exact answer.


Installation
------------

To install the latest release via PyPI using pip::

    pip install interval-binning

The latest development version `can be found on GitHub
<https://github.com/martijnvermaat/interval-binning>`_.


API documentation
-----------------

.. module:: binning

All positions and ranges in this module are zero-based and open-ended,
following standard Python indexing and slicing.

.. autofunction:: assign_bin
.. autofunction:: overlapping_bins
.. autofunction:: containing_bins
.. autofunction:: contained_bins
.. autofunction:: covered_interval


Copyright
---------

This library is licensed under the MIT License, meaning you can do whatever
you want with it as long as all copies include these license terms. The full
license text can be found in the LICENSE.rst file.

See the AUTHORS.txt for for a complete list of copyright holders.
