import numpy as np
from scipy.spatial import Voronoi
from scipy.spatial import distance
from sympy import symbols, Eq, solve
import math
import itertools
import copy
from scipy.spatial import Delaunay
import pandas as pd
import sys
from numpy import linalg as la
from itertools import *
from copy import deepcopy
from collections import defaultdict

from itertools import combinations
from math import prod
from sympy import factorint


        # генерация комбинаций

def generate_factor_combinations(factors):
    if not factors:
        return [[]]
    
    first, *rest = factors
    rest_combinations = generate_factor_combinations(rest)
    
    result = []
    for comb in rest_combinations:
        result.append([first] + comb)
        result.append([first * (comb[0] if comb else 1)] + (comb[1:] if len(comb) > 1 else []))
    
    return result


# находит все возможные разложения числа на простые множители

def compute_factorizations(n):
    
    # находим простые множители числа
    
    factors = [] # список простых множителей
    
    for prime, exp in factorint(n).items():
        factors.extend([prime] * exp)
    
    combinations = []
    combinations = generate_factor_combinations(factors) # генерируем список всех возможных комбинаций
    
    unique = set(tuple(sorted(comb)) for comb in combinations) # убираем дубликаты и сортируем
    unique_combinations = [list(comb) for comb in sorted(unique)]
    
    unique_combinations_4 = []

    for comb in unique_combinations:

        if len(comb) <= 4:
            
            unique_combinations_4.append(comb)

    
    return unique_combinations_4

# добавляем в каждый список 1, чтобы длина каждого списка была = 4
def pad_lists_with_ones(list_of_lists):
    # Проходим по каждому списку в основном списке
    for i in range(len(list_of_lists)):
        # Если длина списка меньше 4
        if len(list_of_lists[i]) < 4:
            # Добавляем в начало столько единиц, сколько не хватает до 4
            list_of_lists[i] = [1] * (4 - len(list_of_lists[i])) + list_of_lists[i]
    return list_of_lists



