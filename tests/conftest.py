"""
Test configuration.
"""


import pytest


@pytest.fixture
def random():
    pass


def pytest_addoption(parser):
    parser.addoption(
        '--repeat-random', metavar='NUM', dest='repeat_random', type=int,
        default=10, help='Number of times to repeat randomized tests')


def pytest_generate_tests(metafunc):
    if 'random' in metafunc.fixturenames:
        repeat_random = metafunc.config.option.repeat_random
        metafunc.fixturenames.append('repeat_count')
        metafunc.parametrize('repeat_count', range(repeat_random))
