"""
Unit tests for the `binning` package using random intervals.
"""


import random

import pytest

import binning


MAX_POSITION = binning.MAX_POSITION
MAX_BIN = binning.MAX_BIN
NUM_INTERVALS = 1000


def random_interval():
    start = random.randint(0, MAX_POSITION)
    stop = random.randint(0, MAX_POSITION + 1)
    if stop < start:
        start, stop = stop, start
    return start, stop


def random_bin():
    return random.randint(0, MAX_BIN)


@pytest.fixture
def interval(random):
    return random_interval()


@pytest.fixture
def intervals(random):
    return set(random_interval() for _ in range(NUM_INTERVALS))


@pytest.fixture
def bin(random):
    return random_bin()


def test_containing(intervals, interval):
    start, stop = interval

    # Intervals completely containing the query interval.
    containing = set((x, y) for x, y in intervals
                     if x <= start and stop <= y)

    # Pre-selection of intervals using binning.
    binned = set((x, y) for x, y in intervals
                 if binning.assign_bin(x, y)
                 in binning.containing_bins(start, stop))

    assert binned.issuperset(containing)


def test_contained(intervals, interval):
    start, stop = interval

    # Intervals completely contained by the query interval.
    contained = set((x, y) for x, y in intervals
                    if start <= x and y <= stop)

    # Pre-selection of intervals using binning.
    binned = set((x, y) for x, y in intervals
                 if binning.assign_bin(x, y)
                 in binning.contained_bins(start, stop))

    assert binned.issuperset(contained)


def test_overlapping(intervals, interval):
    start, stop = interval

    # Intervals overlapping the query interval.
    overlapping = set((x, y) for x, y in intervals
                      if x < stop and start < y)

    # Pre-selection of intervals using binning.
    binned = set((x, y) for x, y in intervals
                 if binning.assign_bin(x, y)
                 in binning.overlapping_bins(start, stop))

    assert binned.issuperset(overlapping)


def test_covered_interval_assign_bin(bin):
    assert binning.assign_bin(*binning.covered_interval(bin)) == bin
