# перебор матриц решетки (float вариант)

import numpy as np
import random
from LLL import *

precision = 3  # точность округления

def normalize_rows(matrix):
    """
    Для каждой строки матрицы:
    - Находит первый ненулевой элемент
    - Если он отрицательный, умножает строку на -1
    """
    matrix = matrix.copy()
    for i in range(matrix.shape[0]):
        row = matrix[i]
        # Находим первый ненулевой элемент
        for val in row:
            if val != 0:
                if val < 0:
                    matrix[i] *= -1
                break  # выходим после первого ненулевого
    return matrix

def canonical_form_by_rows(matrix: np.ndarray) -> np.ndarray:
    """
    Приводит матрицу к каноническому виду через сортировку строк.
    
    matrix — исходная матрица (np.ndarray)
    Возвращает: новая матрица с отсортированными строками
    """

    # Сортируем строки 
    idx = np.lexsort(matrix.T[::-1])  # Инвертируем столбцы для нужного порядка
    return matrix[idx]

def generate_matrix(a = 2):
    # Фиксируем первую строку
    first_row = [1.0, 0.0, 0.0, 0.0]
    
    while True:
        # Генерируем оставшиеся 3 строки
        matrix = [first_row]
        for _ in range(3):
            # Первый элемент — только от -1 до 1
            first = round(random.uniform(0, 1), precision)
            # Остальные — от -a до a
            rest = [round(random.uniform(-a, a), precision) for _ in range(3)]
            row = [first] + rest
            matrix.append(row)
        
        mat = np.array(matrix, dtype=float)
        
        # Проверяем невырожденность
        if not np.isclose(np.linalg.det(mat), 0, atol=1e-8):
            return mat
        
# проверка что длина векторов больше 1

def vector_length(v):
    return np.linalg.norm(v)

def is_valid_matrix(matrix):
    for row in matrix:
        if vector_length(row) < 1:
            return False
    return True



def checkGrid(grid, norm_factor = 2.0, cos_limit = 0.5):
	for i in range(len(grid) - 1):
		vi = grid[i]
		ni = np.linalg.norm(vi)
		for j in range(i+1, len(grid)):
			vj = grid[j]
			nj = np.linalg.norm(vj)

			if ni < 1: return False
			if nj < 1: return False
			if ni / nj > norm_factor: return False
			if nj / ni > norm_factor: return False
   
			cos = vi @ vj / (ni * nj)
			if cos > cos_limit: return False
			if cos < -cos_limit: return False
   
	return True
   

def generateGrids(norm_limit = 2.0, norm_factor = 2.0, cos_limit = 0.5, det_limit = 1.0):
    list_grids = []  # список уникальных базисов
    c = 0  # счётчик хороших матриц
    #k = 0  # общий счётчик
    scale = 10**precision  # Масштаб для округления
    max_attempts = 5000  # Максимальное количество попыток генерации матрицы

    for i in range(max_attempts):
        matrix = generate_matrix(norm_limit)

        if not np.isclose(np.linalg.det(matrix), 0):  # Если матрица невырожденная
        #if is_valid_matrix(matrix):# проверяем что длина векторов больше 1

            # строим LLL базис
            #scale_mat = np.round(matrix * scale).astype(int)
            #grid_int = (scale_mat).astype(int)
            #basis = LLL.reduction(IntegerMatrix.from_matrix(grid_int.tolist()))
            #rows, cols = basis.nrows, basis.ncols
            
            # проверяем уникальность 
            #lll_grid = np.array([[basis[i, j] for j in range(cols)] for i in range(rows)]) 
            #lll_grid = lll_grid * 1/scale # возвращаем float

            lll_grid = LLL(matrix)
            lll_grid = normalize_rows(lll_grid) # делаем все первые элементы неотрицательными
            lll_grid = canonical_form_by_rows(lll_grid) # сортируем строки

            #print(lll_grid, matrix, np.linalg.det(matrix))
            #print(66)
            
            # Проверяем уникальность (без учёта транспонирования)
            #if not any(np.allclose(lll_grid, m) for m in list_grids):
            #    c += 1
            #    list_grids.append(lll_grid)
            
            det = np.linalg.det(lll_grid)
            if det > det_limit and checkGrid(lll_grid, 2.0, 0.2):
                c += 1
                list_grids.append(lll_grid)

        else:
            # Пропускаем эту матрицу
            continue
    return list_grids

def calculateMinMaxDet(list_grids):
    # вычисляем минимальный и максимальный определитель
    # для всех матриц в списке list_grids
    minDet = 1000000000
    maxDet = 0

    for grid in list_grids:
        det = np.abs(np.linalg.det(grid))
        if det < minDet:
            minDet = det
            print(grid, "\nminDet =", minDet)
        if det > maxDet:
            maxDet = det
            print(grid, "\nmaxDet =", maxDet)
    
    return minDet, maxDet