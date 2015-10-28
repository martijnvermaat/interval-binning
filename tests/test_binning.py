"""
Unit tests for the `binning` package.
"""


import pytest

import binning


# Some example intervals with pre-calculated bin numbers.
# http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
INTERVALS = {
    (0, 1): 585,
    (binning.MAX_POSITION, binning.MAX_POSITION + 1): 4680,
    (0, binning.MAX_POSITION + 1): 0,
    (0, 2**17): 585,
    (1, 2**17): 585,
    (0, 2**17 - 1): 585,
    (0, 2**17 + 1): 73,
    (340, 74012): 585,
    (0, 2**18): 73,
    (74012, 173034): 73,
    (423427, 423428): 588,
    (100000, 200000): 73,
    (1000000, 2000000): 9,
    (10000000, 20000000): 1,
    (100000000, 200000000): 0,
    (200000, 1000000): 73,
    (2000000, 10000000): 1,
    (20000000, 100000000): 0,
    (300000000, 300000015): 2873,
    (300000000, 300100015): 359,
    (300000000, 300200015): 359,
    (300000000, 301000015): 44,
    (300000000, 311000015): 5,
    (300000000, 321000015): 5,
    (300000000, 381000015): 0,
    (300000000, 511000015): 0,
    (1200000, 2000000): 74
}


OUT_OF_RANGE_INTERVALS = [
    (-23442, -334),
    (-23442, 334),
    (-23442, 0),
    (-1, -1),
    (-1, 0),
    (5656, binning.MAX_POSITION + 2),
    (-34234, binning.MAX_POSITION + 3432)
]


@pytest.mark.parametrize('start,stop,expected', [
    (start, stop, bin) for (start, stop), bin in INTERVALS.items()])
def test_assign_bin(start, stop, expected):
    assert binning.assign_bin(start, stop) == expected


@pytest.mark.parametrize('start,stop', OUT_OF_RANGE_INTERVALS)
def test_assign__bin_range(start, stop):
    with pytest.raises(binning.OutOfRangeError):
        binning.assign_bin(start, stop)


@pytest.mark.parametrize('start,stop,expected', [
    (0, 1, [(585, 585),
            (73, 73),
            (9, 9),
            (1, 1),
            (0, 0)]),
    (binning.MAX_POSITION, binning.MAX_POSITION + 1, [(4680, 4680),
                                                      (584, 584),
                                                      (72, 72),
                                                      (8, 8),
                                                      (0, 0)]),
    (0, binning.MAX_POSITION + 1, [(585, 4680),
                                   (73, 584),
                                   (9, 72),
                                   (1, 8),
                                   (0, 0)]),
    (0, 2**17 + 1, [(585, 586),
                    (73, 73),
                    (9, 9),
                    (1, 1),
                    (0, 0)]),
    (1200000, 2000000, [(594, 600),
                        (74, 74),
                        (9, 9),
                        (1, 1),
                        (0, 0)]),
    (0, 2**18, [(585, 586),
                (73, 73),
                (9, 9),
                (1, 1),
                (0, 0)]),
    (300000000, 300200015, [(2873, 2875),
                            (359, 359),
                            (44, 44),
                            (5, 5),
                            (0, 0)]),
    (300000000, 301000015, [(2873, 2881),
                            (359, 360),
                            (44, 44),
                            (5, 5),
                            (0, 0)])
])
def test_range_per_level(start, stop, expected):
    assert list(binning.range_per_level(start, stop)) == expected


@pytest.mark.parametrize('start', [0, 1, 2, 2**17 - 2, 2**17 - 1, 2**17,
                                   2**17 + 1, 2**17 + 2, 2**18, 56983456,
                                   binning.MAX_POSITION - 2,
                                   binning.MAX_POSITION - 1,
                                   binning.MAX_POSITION])
def test_range_per_level_zero_length(start):
    assert list(binning.range_per_level(start, start)) == \
        list(binning.range_per_level(start, start + 1))


@pytest.mark.parametrize('start,stop', OUT_OF_RANGE_INTERVALS)
def test_range_per_level_range(start, stop):
    with pytest.raises(binning.OutOfRangeError):
        list(binning.range_per_level(start, stop))


@pytest.mark.parametrize('start,stop,expected', [
    (0, 1,
     [585, 73, 9, 1, 0]),
    (binning.MAX_POSITION, binning.MAX_POSITION + 1,
     [4680, 584, 72, 8, 0]),
    (0, binning.MAX_POSITION + 1,
     list(range(585, 4681)) + list(range(73, 585)) + list(range(9, 73)) +
     list(range(1, 9)) + [0]),
    (0, 2**17 + 1,
     [585, 586, 73, 9, 1, 0]),
    (1200000, 2000000,
     list(range(594, 601)) + [74, 9, 1, 0]),
    (0, 2**18,
     [585, 586, 73, 9, 1, 0]),
    (300000000, 300200015,
     [2873, 2874, 2875, 359, 44, 5, 0]),
    (300000000, 301000015,
     list(range(2873, 2882)) + [359, 360, 44, 5, 0])
])
def test_overlapping_bins(start, stop, expected):
    assert binning.overlapping_bins(start, stop) == expected


@pytest.mark.parametrize('start,stop,expected', [
    (0, 1,
     [585, 73, 9, 1, 0]),
    (binning.MAX_POSITION, binning.MAX_POSITION + 1,
     [4680, 584, 72, 8, 0]),
    (0, binning.MAX_POSITION + 1,
     [0]),
    (0, 2**17 + 1,
     [73, 9, 1, 0]),
    (1200000, 2000000,
     [74, 9, 1, 0]),
    (0, 2**18,
     [73, 9, 1, 0]),
    (300000000, 300200015,
     [359, 44, 5, 0]),
    (300000000, 301000015,
     [44, 5, 0])
])
def test_containing_bins(start, stop, expected):
    assert binning.containing_bins(start, stop) == expected


@pytest.mark.parametrize('start,stop,expected', [
    (0, 1,
     [585]),
    (binning.MAX_POSITION, binning.MAX_POSITION + 1,
     [4680]),
    (0, binning.MAX_POSITION + 1,
     list(range(585, 4681)) + list(range(73, 585)) + list(range(9, 73)) +
     list(range(1, 9)) + [0]),
    (0, 2**17 + 1,
     [585, 586, 73]),
    (1200000, 2000000,
     list(range(594, 601)) + [74]),
    (0, 2**18,
     [585, 586, 73]),
    (300000000, 300200015,
     [2873, 2874, 2875, 359]),
    (300000000, 301000015,
     list(range(2873, 2882)) + [359, 360, 44])
])
def test_contained_bins(start, stop, expected):
    assert binning.contained_bins(start, stop) == expected


@pytest.mark.parametrize('start,stop', INTERVALS.keys())
def test_assign_bin_covered_interval(start, stop):
    bin_start, bin_stop = binning.covered_interval(binning.assign_bin(start,
                                                                      stop))
    assert bin_start <= start and stop <= bin_stop
