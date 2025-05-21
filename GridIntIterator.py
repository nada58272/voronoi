# перебор матриц решетки (целочисленный вариант)

import numpy as np
import itertools
from itertools import combinations
from fpylll import IntegerMatrix, LLL
from tqdm import tqdm
from math import comb

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

# # Список матриц решетки grid
# list_coords = []

# for var in itertools.product(range(-1, 2), repeat=4): 
#     # меньше 3 нельзя - тк при поиске центр региона отбрасываем незаконченные регионы
#     #if var
#     list_coords.append(var)
    
# list_coords.remove((1, 0, 0, 0))
    
# print(len(list_coords), list_coords[12])

# result = []
# for lst in list_coords:
#     first_non_zero = next((x for x in lst if x != 0), None)
#     if first_non_zero is not None and first_non_zero > 0:
#         result.append(lst)

# print(len(result))

# #result = [np.random.rand(4) for _ in range(10)]  # список векторов
# list_grids = []  # список уникальных базисов
# c = 0  # счётчик хороших матриц
# k = 0  # общий счётчик

# # Общее количество комбинаций (для корректного tqdm)
# n = len(result)
# total_combinations = comb(n, 3)

# print(f"Всего комбинаций: {total_combinations}")

# Основной цикл с прогресс-баром
for vectors in tqdm(combinations(result, 3), total=total_combinations, desc="Обработка комбинаций"):
    g = np.vstack(([1, 0, 0, 0], *vectors))  # Собираем векторы в матрицу 4x4
    k += 1
    
    if not np.isclose(np.linalg.det(g), 0):  # Если матрица невырожденная
        # Проверяем уникальность (без учёта транспонирования)
        #k += 1

        #grid_int = (g).astype(int)
        #basis = LLL.reduction(IntegerMatrix.from_matrix(grid_int.tolist()))
        #rows, cols = basis.nrows, basis.ncols
    
        #lll_grid = np.array([[basis[i, j] for j in range(cols)] for i in range(rows)])

        lll_grid = LLL(g)
        
        lll_grid = normalize_rows(lll_grid) # делаем все первые элементы неотрицательными
        lll_grid = canonical_form_by_rows(lll_grid)

        if not any(np.allclose(lll_grid, m) for m in list_grids):
        #if not any(np.allclose(g.T, m) for m in list_grids):
            c += 1
            list_grids.append(lll_grid)

