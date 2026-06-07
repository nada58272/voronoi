"""Тесты разложения определителя на множители (factorization.py)."""

from voronoi4d import compute_factorizations, pad_lists_with_ones


def test_compute_factorizations_simple():
    assert compute_factorizations(12) == [[2, 2, 3], [2, 6], [3, 4], [12]]


def test_compute_factorizations_composite():
    assert compute_factorizations(25) == [[5, 5], [25]]


def test_compute_factorizations_limit_length():
    # длина каждого разложения не превышает 4
    for item in compute_factorizations(32):
        assert len(item) <= 4


def test_pad_lists_with_ones():
    padded = pad_lists_with_ones([[12], [2, 6], [2, 2, 3]])
    assert padded == [
        [1, 1, 1, 12],
        [1, 1, 2, 6],
        [1, 2, 2, 3],
    ]


def test_pad_lists_with_ones_does_not_mutate_input():
    original = [[12]]
    pad_lists_with_ones(original)
    assert original == [[12]]
