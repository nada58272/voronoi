"""Поиск оптимальных подрешёток для периодических раскрасок.

Определитель матрицы перехода M равен индексу подрешётки [Λ:Γ] = |det(M)| —
количеству цветов раскраски (хроматическому числу k).

Основная функция — find_optimal(): для каждого определителя из диапазона
перебирает верхнетреугольные матрицы перехода (эрмитова нормальная форма),
для каждой подрешётки LLL-приводит базис, ищет минимальное нормированное
расстояние d между областями Вороного одного цвета и выбирает матрицу
с максимальным таким расстоянием. Раскраска пригодна при d >= threshold (= 1).
"""

import numpy as np
from itertools import product
from scipy.spatial import distance

from .distances import dist_to_s
from .factorization import compute_factorizations, pad_lists_with_ones
from .io import save_result
from .lll import lll_reduce

# --------------------------------------------------------------------------------


def lattice_points_no_central_symmetry(basis, limits, max_len):
    """Генерирует точки решётки без центрально-симметричных дубликатов.

    Из каждой пары точек (p, -p) остаётся одна: первый ненулевой коэффициент
    должен быть положительным. Если минимальное расстояние между точками
    не превышает max_len, решётка отбрасывается (возвращается только ноль).

    :param basis: базис решётки (numpy.ndarray).
    :param limits: диапазон коэффициентов (-limits, limits).
    :param max_len: диаметр многогранника.
    :return: список точек решётки (включая нулевую).
    """
    points = []

    for coeffs in product(range(-limits, limits + 1), repeat=4):
        if all(c == 0 for c in coeffs):  # пропускаем нулевую точку
            continue

        # каноничность: первый ненулевой коэффициент должен быть положительным
        first_nonzero = next((c for c in coeffs if c != 0), None)
        if first_nonzero is not None and first_nonzero < 0:
            continue

        points.append(np.dot(coeffs, basis))

    points.append(np.array([0, 0, 0, 0]))

    # если минимальное расстояние между точками слишком мало, решётка не подходит
    if len(points) >= 2:
        min_dist = np.min(distance.pdist(points))
        if min_dist <= max_len:
            return [np.array([0, 0, 0, 0])]

    return points


# --------------------------------------------------------------------------------


def find_optimal(det_range, limits, grid, vor4, max_len, precision=3, threshold=1.0,
                 output_file="results.txt", verbose=True):
    """Поиск матриц перехода и расстояний для диапазона определителей.

    :param det_range: range или список определителей (количеств цветов k) для обработки.
    :param limits: границы генерации коэффициентов точек подрешётки.
    :param grid: исходная решётка (numpy.ndarray).
    :param vor4: объект VoronoiPolyhedra после build().
    :param max_len: диаметр центрального многогранника diam(V0).
    :param precision: точность целочисленного масштабирования для LLL.
    :param threshold: минимально допустимое нормированное расстояние d
                      (матрицы с меньшим пропускаются; d >= 1 — пригодная раскраска).
    :param output_file: файл для записи результатов.
    :param verbose: печатать прогресс.
    :return: словари det_dist, det_center, det_mat (ключ — определитель).
    """
    centers_dist = {}  # кэш: ключ — координаты точки, значение — расстояние

    det_dist = {}  # ключ — определитель, значение — расстояние
    det_center = {}  # ключ — определитель, значение — координаты точки s
    det_mat = {}  # ключ — определитель, значение — матрица перехода

    for det in det_range:
        if verbose:
            print("\r                                                          ")
            print("---------------------------------")
            print("det:", det)

        # все варианты диагоналей матрицы перехода с данным определителем
        list_diag_el = pad_lists_with_ones(compute_factorizations(det))

        # среди расстояний ищем максимальное по всем матрицам mat
        mat_dist = {}
        mat_center = {}
        list_mats = []
        index = 0

        for diag_el in list_diag_el:
            mat = np.diag(np.array(diag_el, dtype=float))
            max_num_col1, max_num_col2, max_num_col3 = diag_el[1], diag_el[2], diag_el[3]

            num_iterations = max_num_col1 * max_num_col2 ** 2 * max_num_col3 ** 3
            iteration = 0

            if verbose:
                print("\r                                                          ")
                print("▶ diag factors:", *diag_el, "   iters:", num_iterations)

            # перебираем наддиагональные элементы (эрмитова нормальная форма)
            for indices in product(range(max_num_col3), range(max_num_col3), range(max_num_col2),
                                   range(max_num_col3), range(max_num_col2), range(max_num_col1)):

                mat[2][3] = indices[0]
                mat[1][3] = indices[1]
                mat[1][2] = indices[2]
                mat[0][3] = indices[3]
                mat[0][2] = indices[4]
                mat[0][1] = indices[5]

                iteration += 1
                if verbose and iteration % 500 == 0:
                    print("\r[", int(10000 * iteration / num_iterations) / 100, "% ]",
                          "   iter:", iteration, end="")

                # базис подрешётки и его LLL-приведение
                sub_grid = np.dot(mat, grid)
                sub_grid_lll = lll_reduce(sub_grid, precision=precision)

                centers = lattice_points_no_central_symmetry(sub_grid_lll, limits, vor4.max_len)

                # решётка отброшена (слишком близкие точки)
                if len(centers) == 1 and np.array_equal(centers[0], np.array([0, 0, 0, 0])):
                    continue

                # минимальное расстояние для данной матрицы mat
                min_dist_mat = float("inf")
                min_center = None

                for center in centers:
                    center_key = tuple(center)

                    if center_key in centers_dist:
                        dist = centers_dist[center_key]
                    else:
                        s = 0.5 * center
                        dist = dist_to_s(vor4, s, max_len)
                        centers_dist[center_key] = dist

                    if dist < min_dist_mat:
                        min_dist_mat = dist
                        min_center = center.copy()

                    if min_dist_mat < threshold:
                        break

                # расстояние меньше порога — матрица не подходит
                if min_dist_mat < threshold:
                    continue

                if verbose:
                    print("\r", mat, "                      \n", min_dist_mat, min_center)

                # сохраняем значения для mat (копию: mat мутируется в цикле!)
                list_mats.append(mat.copy())
                mat_dist[index] = min_dist_mat
                mat_center[index] = min_center
                index += 1

        # максимальное среди минимальных расстояний для данного определителя
        if mat_dist:
            best_index, best_dist = max(mat_dist.items(), key=lambda item: item[1])
            det_dist[det] = best_dist
            det_center[det] = mat_center[best_index]
            det_mat[det] = list_mats[best_index]

            save_result(grid, det, list_mats[best_index], mat_center[best_index], best_dist,
                        output_file=output_file)

    return det_dist, det_center, det_mat
