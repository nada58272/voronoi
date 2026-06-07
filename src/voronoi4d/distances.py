"""Расчёт расстояний от точек до центрального многогранника Вороного.

Основная функция — dist_to_s(): расстояние от точки s до центрального
многогранника через каскад проекций (3D-грань → 2D-грань → ребро → вершина).
"""

import math
import warnings

import numpy as np
from scipy.spatial import distance

CHECK_DIST = True  # проверка расстояний по теореме Пифагора (можно отключить: установить False)
TOL_SIMPLEX = 1e-9  # допуск поиска симплекса в триангуляции
TOL_NEAREST = 1e-7  # допуск сравнения расстояний до ближайших вершин

# --------------------------------------------------------------------------------


def check_dist(dist1, dist2):
    """Сверяет два квадрата расстояний (контроль по теореме Пифагора).

    При расхождении выдаёт предупреждение через модуль warnings (не прерывает счёт).
    """
    if not math.isclose(dist1, dist2, abs_tol=1e-9):
        warnings.warn(
            f"расхождение квадратов расстояний (теорема Пифагора): {dist1 - dist2}",
            stacklevel=2,
        )


# --------------------------------------------------------------------------------


def find_faces_from_nearest_vertices(s, central, vertex_to_faces):
    """Находит грани, связанные с ближайшими к точке s вершинами.

    :param s: координаты точки (np.array).
    :param central: список вершин центрального многогранника.
    :param vertex_to_faces: для каждой вершины — список индексов содержащих её граней.
    :return: множество индексов граней, общих для всех ближайших вершин.
    """
    # два прохода, чтобы не терять со-ближайшую вершину при численном равенстве:
    # сначала находим минимальное расстояние, затем собираем все вершины в его пределах
    dists = [distance.euclidean(s, central[index]) for index in range(len(central))]
    min_dist_vert_to_s = min(dists)  # минимальное расстояние до точки s

    # индексы вершин, со-ближайших к s в пределах допуска
    min_dist_to_s_list = [
        index
        for index, dist_vert_to_s in enumerate(dists)
        if np.allclose(dist_vert_to_s - min_dist_vert_to_s, 0, atol=TOL_NEAREST)
    ]

    # пересечение множеств граней, связанных с ближайшими вершинами
    selected_lists = [set(vertex_to_faces[i]) for i in min_dist_to_s_list]
    return set.intersection(*selected_lists)


# --------------------------------------------------------------------------------


def dist_to_s(vor4, s, max_len):
    """Нормированное расстояние от точки s до центрального многогранника V0.

    s — середина отрезка между началом координат и центром соседней области
    Вороного подрешётки. По лемме Иванова минимальное расстояние D между
    соседними областями реализуется в точках, симметричных относительно s,
    поэтому D равно удвоенному расстоянию от s до V0. Возвращается
    нормированное запрещённое расстояние d = 2*D / diam(V0) = dist * 2 / max_len.

    Проекция ищется каскадом: на 3-мерную грань, затем (если проекция вне
    многогранника) на её 2-мерные грани, затем на рёбра и вершины.

    Ранний выход: как только нормированное расстояние очередной грани оказывается
    меньше 1, функция сразу возвращает ТЕКУЩИЙ минимум (не обязательно глобальный
    по всем граням). Для алгоритма раскраски этого достаточно — важен лишь факт
    d < 1: раскраска с таким запрещённым расстоянием заведомо непригодна.

    :param vor4: объект VoronoiPolyhedra после build().
    :param s: координаты точки (np.array).
    :param max_len: диаметр центрального многогранника diam(V0).
    :return: нормированное расстояние d.
    """
    polyhedrons = vor4.polyhedrons

    min_dist_to_pol = float("inf")  # минимальное расстояние до центрального многогранника

    # грани, которым принадлежит ближайшая к s вершина
    nearest_faces = find_faces_from_nearest_vertices(s, vor4.central, vor4.vertex_to_faces)

    def update_min_distance():
        nonlocal min_dist_to_pol

        if dist < min_dist_to_pol:
            min_dist_to_pol = dist

    # рассматриваем только грани, которым принадлежит ближайшая к s вершина
    for i in nearest_faces:
        # проекция на 3-мерную грань
        d0 = polyhedrons[i].normal @ (s - polyhedrons[i].center)
        coord0 = s - d0 * polyhedrons[i].normal
        simplex = vor4.delaunay.find_simplex(coord0, tol=TOL_SIMPLEX)

        d0_squared = d0 * d0

        if simplex != -1:  # проекция принадлежит центральному многограннику
            dist = abs(d0)
            update_min_distance()
            continue

        for face2d in polyhedrons[i].faces:
            # проекция на 2-мерную грань
            d1 = face2d.normal @ (coord0 - face2d.center)
            coord1 = coord0 - d1 * face2d.normal
            simplex = vor4.delaunay.find_simplex(coord1, tol=TOL_SIMPLEX)

            d1_squared = d1 * d1

            if simplex != -1:  # проекция принадлежит центральному многограннику
                dist = distance.euclidean(s, coord1)

                if CHECK_DIST:
                    check_dist(dist * dist, d0_squared + d1_squared)

                update_min_distance()
                continue

            for edge in face2d.edges:
                # проекция на ребро
                d2 = edge.normal @ (coord1 - edge.center)
                coord2 = coord1 - d2 * edge.normal
                simplex = vor4.delaunay.find_simplex(coord2, tol=TOL_SIMPLEX)

                d2_squared = d2 * d2

                if simplex != -1:  # проекция принадлежит центральному многограннику
                    dist = distance.euclidean(s, coord2)

                    if CHECK_DIST:
                        check_dist(dist * dist, d0_squared + d1_squared + d2_squared)

                    update_min_distance()
                else:
                    # проекция вне ребра — берём ближайшую вершину ребра
                    d3 = distance.euclidean(coord2, edge.vertex1)
                    d4 = distance.euclidean(coord2, edge.vertex2)

                    if d3 < d4:
                        dist = distance.euclidean(s, edge.vertex1)
                        d34_squared = d3 * d3
                    else:
                        dist = distance.euclidean(s, edge.vertex2)
                        d34_squared = d4 * d4

                    if CHECK_DIST:
                        check_dist(dist * dist, d0_squared + d1_squared + d2_squared + d34_squared)

                    update_min_distance()

        # если нормированное расстояние уже < 1, дальше можно не считать
        if min_dist_to_pol * 2 / max_len < 1:
            return min_dist_to_pol * 2 / max_len

    return min_dist_to_pol * 2 / max_len
