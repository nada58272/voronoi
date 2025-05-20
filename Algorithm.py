import numpy as np
from itertools import product
#from fpylll import IntegerMatrix, LLL

from Distances import *
from Factorization import *
from LLL import *

# Файл для логирования результатов
output_file = "results.txt"

def write_result_to_file(grid, det, mat, center, distance):
    """
    Записывает результаты в файл.
    :param grid: исходная решётка (numpy.ndarray)
    :param det: определитель матрицы
    :param mat: текущая матрица перехода (numpy.ndarray)
    :param center: координаты точки (numpy.ndarray)
    :param distance: расстояние до многогранника
    """
    with open(output_file, 'a') as f:  # 'a' - append mode
        f.write(f"Grid: {grid.tolist()}\n")
        f.write(f"Determinant: {det}\n")
        f.write(f"Matrix:\n{mat.tolist()}\n")
        f.write(f"Center: {center.tolist()}\n")
        f.write(f"Distance: {distance:.6f}\n")
        f.write("-" * 50 + "\n")

#--------------------------------------------------------------------------------

def lattice_points_no_central_symmetry(basis, limits, coeff_dist, max_len):
    """
    Генерирует точки решетки, заданной базисом, без центральносимметричных точек,
    сразу исключая симметричные коэффициенты.
    
    :param limits: Целое число, определяющее диапазон коэффициентов для генерации точек (-limits, limits).
    :param coeff_dist: коэффициент, количество диаметров многогранника для расчета r_max
    :r_max: максимальное расстояние, дальше которого нет смысла рассматривать()
    :return: Список коэффициентов без центральносимметричных точек.  
    """
    
    # Инициализируем список для хранения точек
    points = []
    r_max = coeff_dist * max_len
    

    # Генерируем все возможные целочисленные комбинации коэффициентов
    for coeffs in product(range(-limits, limits + 1), repeat=4):
        if all(c == 0 for c in coeffs):  # Пропускаем нулевую точку
            continue
        
        # Проверяем "каноничность" коэффициентов
        # Каноническая форма: первый ненулевой коэффициент должен быть положительным
        first_nonzero_index = next((i for i, c in enumerate(coeffs) if c != 0), None)
        
        if first_nonzero_index is not None and coeffs[first_nonzero_index] < 0:
            continue  # Пропускаем, если первый ненулевой коэффициент отрицательный

        # Вычисляем точку решетки
        point = np.dot(coeffs, basis)
    
        # Предварительная фильтрация: ограничение по максимальной координате
        max_coord = int(np.ceil(r_max))  # Максимально возможное значение координаты
        
        # Проверяем, что координаты точки находятся в допустимых пределах
        if any(abs(c) >= max_coord for c in point):  # Отсеиваем точки, которые заведомо далеки
            continue
        
        points.append(point)

    # Минимальное расстояние между точками
    if len(points) >= 2:
        dist_matrix = distance.pdist(points)
        min_dist = np.min(dist_matrix)

        if min_dist <= max_len:
            points = []
            points.append(np.array([0, 0, 0, 0]))
            return points
    
    return points

#--------------------------------------------------------------------------------

def find_optimal(det_range, limits, grid, vor4, max_len):
    """
    Основная функция для выполнения алгоритма поиска матриц и расстояний.

    :param det_range: range или список определителей для обработки.
    :param limits: целое число, задающее границы генерации коэффициентов.
    :param grid: исходная решётка (numpy.ndarray).
    :param vor4: объект с многогранниками Вороного.
    :return: словари det_dist, det_center, det_mat с результатами.
    """
    
    # Инициализация словарей для хранения результатов
    centers_dist = dict() 

    det_dist = dict() # словарь где ключ - определитель, значение - расстояние
    det_center = dict() # словарь где ключ - определитель, значение - координаты точки s
    det_mat = dict()

    for det in det_range:

        print("\r                                                          ")
        print("---------------------------------")
        print("det:", det)
        
        list_all_factorizations = compute_factorizations(det)
        list_diag_el = pad_lists_with_ones(list_all_factorizations)
        
        # среди этих расстояний надо найти максимальное для каждой матрицы mat
        mat_dist = dict()
        mat_center = dict()
        list_mats = []
        index = 0
        
        for diag_el in list_diag_el:
            mat = np.array([[diag_el[0], 0, 0, 0], 
                            [0, diag_el[1], 0, 0], 
                            [0, 0, diag_el[2], 0], 
                            [0, 0, 0, diag_el[3]]], float)
            max_num_col3 = diag_el[3]
            max_num_col2 = diag_el[2]
            max_num_col1 = diag_el[1]
            
            num_iterations = diag_el[1] * diag_el[2] * diag_el[2] * diag_el[3] * diag_el[3] * diag_el[3]
            iter = 0
            
            print("\r                                                          ")
            print("▶ diag factors:", diag_el[1], diag_el[2], diag_el[3], "   iters:", num_iterations)

            for indices in product(range(max_num_col3), range(max_num_col3), range(max_num_col2),
                                   range(max_num_col3), range(max_num_col2), range(max_num_col1)):

                mat[2][3] = indices[0]
                mat[1][3] = indices[1]
                mat[1][2] = indices[2]
                mat[0][3] = indices[3]
                mat[0][2] = indices[4]
                mat[0][1] = indices[5]

                iter += 1
                if iter % 500 == 0:
                    print("\r[", int(10000 * iter / num_iterations) / 100, "% ]", "   iter:", iter, end = '')
                    
                sub_grid = np.dot(mat, grid)

                #sub_grid_int = (sub_grid).astype(int)
                #basis = LLL.reduction(IntegerMatrix.from_matrix(sub_grid_int.tolist()))
                #rows, cols = basis.nrows, basis.ncols
                #sub_grid_LLL = np.array([[basis[i, j] for j in range(cols)] for i in range(rows)])
                
                sub_grid_LLL = LLL(sub_grid, 0.75)
               # print(sub_grid_LLL)
                centers = lattice_points_no_central_symmetry(sub_grid_LLL, limits, 3, vor4.max_len)
                min_dist_mat = 2 * vor4.max_len

                if len(centers) == 1 and np.array_equal(centers[0], np.array([0, 0, 0, 0])):
                    continue

                for center in centers:
                    center_key = tuple(center)

                    if center_key in centers_dist:
                        dist = centers_dist[center_key]
                    else:
                        s = 0.5 * center
                        dist, coords, ind = dist_to_s(vor4.polyhedrons, s, vor4, max_len)
                        # добавляем новое значение в centers_dist
                        centers_dist[center_key] = dist

                    if dist < min_dist_mat:
                        min_dist_mat = dist
                        min_center = center.copy()
                        min_mat = mat.copy()
                        min_sub_grid_LLL = sub_grid_LLL.copy()
                        min_sub_grid = sub_grid.copy()

                    if dist > min_dist_mat + vor4.max_len: continue

                    if min_dist_mat < 1: break

                if min_dist_mat < 1: continue
                
                print("\r", min_dist_mat, min_center, min_mat)#, min_sub_grid_LLL, min_sub_grid)#, s, centers)

                # сохраняем значения для mat
                list_mats.append(mat)

                mat_dist[index] = min_dist_mat
                mat_center[index] = min_center
                index += 1
                                    
        # находим максимальный среди минимальных элементов по расстоянию для данного определителя
        if mat_dist:
            max_pair = max(mat_dist.items(), key=lambda item: item[1])
            det_dist[det] = max_pair[1] # ключ - определитель, значение - расстояние
            det_center[det] = mat_center[max_pair[0]] # ключ - определитель, значение - координаты центра
            det_mat[det] = list_mats[max_pair[0]] # ключ - определитель, значение - матрица mat
        
            write_result_to_file(grid, det, list_mats[max_pair[0]], mat_center[max_pair[0]], max_pair[1])
    return det_dist, det_center, det_mat

