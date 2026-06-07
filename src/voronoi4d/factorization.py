"""Разложение определителя на множители.

Для заданного определителя строятся все варианты разложения на ≤4 множителя —
они становятся диагональными элементами матриц перехода к подрешёткам.
"""

from sympy import factorint

# --------------------------------------------------------------------------------


def generate_factor_combinations(factors):
    """Генерирует все комбинации группировки списка простых множителей.

    :param factors: список простых множителей (с повторениями).
    :return: список комбинаций (каждая — список множителей).
    """
    if not factors:
        return [[]]

    first, *rest = factors
    rest_combinations = generate_factor_combinations(rest)

    result = []
    for comb in rest_combinations:
        result.append([first] + comb)
        result.append([first * (comb[0] if comb else 1)] + (comb[1:] if len(comb) > 1 else []))

    return result


def compute_factorizations(n):
    """Находит все разложения числа n на множители длиной не более 4.

    :param n: целое число (определитель).
    :return: отсортированный список уникальных разложений.
    """
    # простые множители числа (с учётом кратности)
    factors = []
    for prime, exp in factorint(n).items():
        factors.extend([prime] * exp)

    combinations = generate_factor_combinations(factors)

    # убираем дубликаты и сортируем
    unique = set(tuple(sorted(comb)) for comb in combinations)
    unique_combinations = [list(comb) for comb in sorted(unique)]

    # оставляем только разложения длиной не более 4
    return [comb for comb in unique_combinations if len(comb) <= 4]


def pad_lists_with_ones(list_of_lists):
    """Дополняет каждый список единицами в начале до длины 4.

    :param list_of_lists: список разложений.
    :return: список разложений, каждое длиной ровно 4.
    """
    return [[1] * (4 - len(lst)) + lst if len(lst) < 4 else lst for lst in list_of_lists]
