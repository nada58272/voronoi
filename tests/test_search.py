"""Интеграционные тесты поиска подрешёток (search.py)."""

import numpy as np
import pytest

from voronoi4d import (
    VoronoiPolyhedra,
    dist_to_s,
    find_optimal,
    lattice_points_no_central_symmetry,
    lll_reduce,
)

# решётка D4 (шахматная): её 24-ячейка компактна относительно подрешёток
D4_GRID = np.array([
    [1, 1, 0, 0],
    [1, -1, 0, 0],
    [0, 1, -1, 0],
    [0, 0, 1, -1],
], dtype=float)


@pytest.fixture(scope="module")
def vor_d4():
    vor = VoronoiPolyhedra(D4_GRID)
    vor.build(verbose=False)
    return vor


def test_find_optimal_rejects_too_dense_sublattices(vor, tmp_path):
    """Для Z^4 и det=4 минимальный вектор подрешётки короче диаметра — всё отбрасывается."""
    grid = np.eye(4, dtype=float)
    output_file = str(tmp_path / "results.txt")

    det_dist, det_center, det_mat = find_optimal(
        range(4, 5), 1, grid, vor, vor.max_len,
        threshold=0.0, output_file=output_file, verbose=False,
    )

    assert det_dist == {}
    assert det_center == {}
    assert det_mat == {}


def test_find_optimal_smoke_d4(vor_d4, tmp_path):
    """Малый end-to-end прогон: det=16 для решётки D4.

    Заодно регрессия: сохранённая матрица должна воспроизводить сохранённое
    расстояние. Раньше в list_mats сохранялась ссылка на мутируемую матрицу,
    и все сохранённые матрицы перетирались её последним состоянием.
    """
    output_file = str(tmp_path / "results.txt")

    det_dist, det_center, det_mat = find_optimal(
        range(16, 17), 1, D4_GRID, vor_d4, vor_d4.max_len,
        threshold=0.0, output_file=output_file, verbose=False,
    )

    assert 16 in det_dist and 16 in det_center and 16 in det_mat

    # матрица перехода верхнетреугольная с определителем 16
    mat = det_mat[16]
    assert np.allclose(mat, np.triu(mat))
    assert np.isclose(np.linalg.det(mat), 16.0)

    # результат записан в файл
    content = open(output_file).read()
    assert "Determinant: 16" in content

    # регрессия mat.copy(): пересчитываем расстояние по сохранённой матрице
    sub_grid = np.dot(mat, D4_GRID)
    sub_grid_lll = lll_reduce(sub_grid, precision=3)
    centers = lattice_points_no_central_symmetry(sub_grid_lll, 1, vor_d4.max_len)

    recomputed = min(
        dist_to_s(vor_d4, 0.5 * center, vor_d4.max_len) for center in centers
    )

    assert np.isclose(recomputed, det_dist[16], atol=1e-9), \
        "Сохранённая матрица не соответствует сохранённому расстоянию"
