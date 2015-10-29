"""
A Python implementation of the interval binning scheme.

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

All positions and ranges in this module are zero-based and open-ended,
following standard Python indexing and slicing.

.. moduleauthor:: Martijn Vermaat <martijn@vermaat.name>

.. Licensed under the MIT license, see the LICENSE.rst file.
"""


# TODO: Implement the extended binning scheme (for positions > 2^29).
# TODO: Be more flexible in the binning scheme to use.


from __future__ import unicode_literals

from itertools import dropwhile


__version_info__ = ('1', '0', '1', 'dev')
__date__ = '28 Oct 2015'


__version__ = '.'.join(__version_info__)
__author__ = 'Martijn Vermaat'
__contact__ = 'martijn@vermaat.name'
__homepage__ = 'https://github.com/martijnvermaat/binning'


# Standard scheme used by the UCSC Genome Browser.
BIN_OFFSETS = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0]
SHIFT_FIRST = 17
SHIFT_NEXT = 3
MAX_POSITION = pow(2, 29) - 1
MAX_BIN = BIN_OFFSETS[0] + (MAX_POSITION >> SHIFT_FIRST)


class OutOfRangeError(Exception):
    """
    Exception that is raised on seeing an invalid bin number or a position or
    interval exceeding the range of the binning scheme.
    """
    pass


def range_per_level(start, stop):
    """
    Given an interval `start:stop`, make an iterator that returns for each
    level the first and last bin overlapping the interval, starting with the
    smallest bins.

    Zero-length intervals are binned according to the one-position interval
    following it.

    Algorithm by `Jim Kent
    <http://genomewiki.ucsc.edu/index.php/Bin_indexing_system>`_.

    :arg int start, stop: Interval positions (zero-based, open-ended).

    :return: Iterator yielding tuples `first, last` being the first and last
      bin overlapping with `start:stop`. The tuples are ordered according to
      the bin size of the levels, starting with the smallest bins.
    :rtype: iterator(tuple(int, int))

    :raise OutOfRangeError: If `start:stop` exceeds the range of the binning
      scheme.
    """
    if start < 0 or stop > MAX_POSITION + 1:
        raise OutOfRangeError(
            'Interval %d:%d is out of range (maximum position is %d)'
            % (start, stop, MAX_POSITION))

    # Note that we treat the zero-length interval `x:x` as `x:x+1`.
    start_bin = start
    stop_bin = max(start, stop - 1)

    start_bin >>= SHIFT_FIRST
    stop_bin >>= SHIFT_FIRST

    for offset in BIN_OFFSETS:
        yield offset + start_bin, offset + stop_bin
        start_bin >>= SHIFT_NEXT
        stop_bin >>= SHIFT_NEXT


def assign_bin(start, stop):
    """
    Given an interval `start:stop`, return the smallest bin in which it fits.

    :arg int start, stop: Interval positions (zero-based, open-ended).

    :return: Smallest bin containing `start:stop`.
    :rtype: int

    :raise OutOfRangeError: If `start:stop` exceeds the range of the binning
      scheme.
    """
    try:
        return next(dropwhile(lambda r: r[0] != r[1],
                              range_per_level(start, stop)))[0]
    except StopIteration:
        raise Exception('An unexpected error occured in assigning a bin')


def overlapping_bins(start, stop=None):
    """
    Given an interval `start:stop`, return bins for intervals *overlapping*
    `start:stop` by at least one position. The order is according to the bin
    level (starting with the smallest bins), and within a level according to
    the bin number (ascending).

    :arg int start, stop: Interval positions (zero-based, open-ended). If
      `stop` is not provided, the interval is assumed to be of length 1
      (equivalent to `stop = start + 1`).

    :return: All bins for intervals overlapping `start:stop`, ordered first
      according to bin level (ascending) and then according to bin number
      (ascending).
    :rtype: list(int)

    :raise OutOfRangeError: If `start:stop` exceeds the range of the binning
      scheme.
    """
    if stop is None:
        stop = start + 1

    return [bin
            for first, last in range_per_level(start, stop)
            for bin in range(first, last + 1)]


def containing_bins(start, stop=None):
    """
    Given an interval `start:stop`, return bins for intervals completely
    *containing* `start:stop`. The order is according to the bin level
    (starting with the smallest bins), and within a level according to the bin
    number (ascending).

    :arg int start, stop: Interval positions (zero-based, open-ended). If
      `stop` is not provided, the interval is assumed to be of length 1
      (equivalent to `stop = start + 1`).

    :return: All bins for intervals containing `start:stop`, ordered first
      according to bin level (ascending) and then according to bin number
      (ascending).
    :rtype: list(int)

    :raise OutOfRangeError: If `start:stop` exceeds the range of the binning
      scheme.
    """
    if stop is None:
        stop = start + 1

    max_bin = assign_bin(start, stop)
    return [bin for bin in overlapping_bins(start, stop) if bin <= max_bin]


def contained_bins(start, stop=None):
    """
    Given an interval `start:stop`, return bins for intervals completely
    *contained by* `start:stop`. The order is according to the bin level
    (starting with the smallest bins), and within a level according to the bin
    number (ascending).

    :arg int start, stop: Interval positions (zero-based, open-ended). If
      `stop` is not provided, the interval is assumed to be of length 1
      (equivalent to `stop = start + 1`).

    :return: All bins for intervals contained by `start:stop`, ordered first
      according to bin level (ascending) and then according to bin number
      (ascending).
    :rtype: list(int)

    :raise OutOfRangeError: If `start:stop` exceeds the range of the binning
      scheme.
    """
    if stop is None:
        stop = start + 1

    min_bin = assign_bin(start, stop)
    return [bin for bin in overlapping_bins(start, stop) if bin >= min_bin]


def covered_interval(bin):
    """
    Given a bin number `bin`, return the interval covered by this bin.

    :arg int bin: Bin number.

    :return: Tuple of `start, stop` being the zero-based, open-ended interval
      covered by `bin`.
    :rtype: tuple(int)

    :raise OutOfRangeError: If bin number `bin` exceeds the maximum bin
      number.
    """
    if bin < 0 or bin > MAX_BIN:
        raise OutOfRangeError(
            'Invalid bin number %d (maximum bin number is %d)'
            % (bin, MAX_BIN))

    shift = SHIFT_FIRST
    for offset in BIN_OFFSETS:
        if offset <= bin:
            return bin - offset << shift, bin + 1 - offset << shift
        shift += SHIFT_NEXT
