import unittest
import numpy as np
import os
import sys

# Попытка импорта fpylll
try:
    import fpylll
    FPYLLL_AVAILABLE = True
except ImportError:
    FPYLLL_AVAILABLE = False
    print("WARNING: Библиотека 'fpylll' не найдена. Тесты для Algorithm.py будут пропущены.\n")

# Импорт модулей проекта

try:
    from Factorization import compute_factorizations, pad_lists_with_ones, generate_factor_combinations
    from Polyhedra import VoronoiPolyhedra, Face2D, Edge2D, Face3D
    from Distances import dist_to_s, find_faces_from_nearest_vertices, checkDist
    from Algorithm import lattice_points_no_central_symmetry
except ImportError as e:
    print(f"Ошибка импорта: {e}")
    print("Убедитесь, что файлы Factorization.py, Polyhedra.py, Distances.py, Algorithm.py находятся в одной папке.")
    exit(1)

    

class TestFactorization(unittest.TestCase):
    """
    Тестирование модуля разложения чисел (Factorization.py).
    """
    
    # проверка разложения числа на множители
    def test_compute_factorizations_simple(self):
        res = compute_factorizations(12)
        self.assertEqual(res, [[2, 2, 3], [2, 6], [3, 4], [12]])
    
    # проверка разложения числа на множители
    def test_compute_factorizations_composite(self):
        res = compute_factorizations(25)
        self.assertEqual(res, [[5, 5], [25]])

    # проверка длины разложения чисел
    def test_compute_factorizations_limit_length(self):
        res = compute_factorizations(32)
        for item in res:
            self.assertLessEqual(len(item), 4)
            
    # проверка дополения единицами до 4 элементов
    def test_pad_lists_with_ones(self):
        input_list = [[12], [2, 6], [2, 2, 3]]
        padded = pad_lists_with_ones(input_list)
        expected = [
            [1, 1, 1, 12],
            [1, 1, 2, 6],
            [1, 2, 2, 3]
        ]
        self.assertEqual(padded, expected)

        

class TestPolyhedra(unittest.TestCase):
    """
    Тестирование построения геометрии Вороного (Polyhedra.py).
    """
    
    @classmethod
    def setUpClass(cls):
        print("\n[Polyhedra] Построение диаграммы Вороного для тестов...")
        # Используем единичную матрицу Z^4
        cls.grid = np.eye(4, dtype=float)
        cls.vor = VoronoiPolyhedra(cls.grid)
        cls.vor.build()

    def test_voronoi_build_integrity(self):
        self.assertTrue(len(self.vor.coords4) > 0, "Координаты не сгенерированы")
        self.assertIsNotNone(self.vor.central_region_index, "Центральный регион не найден")
        self.assertTrue(len(self.vor.faces_3d) > 0, "3D грани не найдены")
        self.assertTrue(len(self.vor.polyhedrons) > 0, "Объекты полиэдров не созданы")

    def test_polyhedron_geometry(self):
        poly = self.vor.polyhedrons[0]
        self.assertEqual(len(poly.normal), 4)
        self.assertAlmostEqual(np.linalg.norm(poly.normal), 1.0, places=5, msg="Нормаль не единичная")
        self.assertTrue(len(poly.faces) > 0)
        self.assertIsInstance(poly, Face3D) 

    def test_max_len_calculation(self):
        self.assertGreater(self.vor.max_len, 0)

    def test_triangulation(self):
        self.assertIsNotNone(self.vor.delaunay)

        
    def test_faces_and_polyhedrons(self):
        # Проверяем, что грани были найдены
        self.assertTrue(len(self.vor.faces_3d) > 0)
        self.assertTrue(len(self.vor.polyhedrons) > 0)
        
        # Для 4D гиперкуба должно быть 8 кубических граней (ячеек)
        self.assertEqual(len(self.vor.polyhedrons), 8) 

        
        
    def test_central_region_found(self):
        self.assertIsNotNone(self.vor.central_region_index)
        self.assertTrue(len(self.vor.central) > 0)
        # Для кубической решетки центральная точка (0,0,0,0) должна быть в coords4
        origin = [0.0, 0.0, 0.0, 0.0]
        self.assertIn(origin, self.vor.coords4)

       
        
        

class TestDistances(unittest.TestCase):
    """
    Тестирование расчета расстояний (Distances.py).
    """

    @classmethod
    def setUpClass(cls):
        cls.grid = np.eye(4, dtype=float)
        cls.vor = VoronoiPolyhedra(cls.grid)
        cls.vor.build()

    def test_lattice_points_no_central_symmetry(self):
        points = lattice_points_no_central_symmetry(np.eye(4), 1, 2, self.vor.max_len)
        
        self.assertIsInstance(points, list)
        has_zero = any(np.array_equal(p, np.zeros(4)) for p in points)
        self.assertTrue(has_zero, "Нулевая точка должна присутствовать в списке")

        points_list = [list(p) for p in points]
        if [1.0, 0.0, 0.0, 0.0] in points_list:
             self.assertNotIn([-1.0, 0.0, 0.0, 0.0], points_list)


    def test_checkDist(self):
        try:
            checkDist(1.0, 1.0000000001)
            checkDist(1.0, 2.0)
        except Exception as e:
            self.fail(f"checkDist упала с ошибкой: {e}")
            
            
            
            
    def test_dist_orthogonal(self):
        """
        Тест 1: Точка строго напротив центра грани.
        Случай: Точка (1.0, 0, 0, 0).
        Ближайшая грань куба находится при x = 0.5.
        Ожидаемое расстояние: 1.0 - 0.5 = 0.5.
        """
        s = np.array([1.0, 0.0, 0.0, 0.0])
        
        # max_len для куба = 2.0. Функция возвращает: raw_dist * 2 / max_len
        # Значит ожидаем: 0.5 * 1.0 = 0.5
        result = dist_to_s(self.vor, s, self.vor.max_len)
        
        self.assertAlmostEqual(result, 0.5, places=5, 
                               msg="Ошибка при расчете расстояния по прямой (ортогонально)")

    def test_dist_diagonal_vertex(self):
        """
        Тест 2: Точка на диагонали, за вершиной куба.
        Вершина куба имеет координаты (0.5, 0.5, 0.5, 0.5).
        Берем точку s = (1.0, 1.0, 1.0, 1.0).
        Расстояние должно считаться от s до Вершины.
        """
        vertex = np.array([0.5, 0.5, 0.5, 0.5])
        s = np.array([1.0, 1.0, 1.0, 1.0])
        
        # Ручной расчет евклидова расстояния
        expected_raw_dist = np.linalg.norm(s - vertex) # sqrt(0.5^2 * 4) = sqrt(1) = 1.0
        
        result = dist_to_s(self.vor, s, self.vor.max_len)
        
        self.assertAlmostEqual(result, expected_raw_dist, places=5, 
                               msg="Ошибка при расчете расстояния до вершины (диагональ)")

    def test_point_on_surface(self):
        """
        Тест 3: Точка лежит прямо на грани.
        Расстояние должно быть 0.
        """
        s = np.array([0.5, 0.0, 0.0, 0.0]) # Эта точка лежит на грани x=0.5
        
        result = dist_to_s(self.vor, s, self.vor.max_len)
        
        self.assertAlmostEqual(result, 0.0, places=5, 
                               msg="Расстояние для точки на поверхности должно быть 0")

    def test_point_slightly_offset(self):
        """
        Тест 4: Проверка проекции. Смещаем точку и по нормали, и вдоль грани.
        Точка (1.7, 0.1, 0.1, 0.1).
        Проекция на куб должна упасть на грань x=0.5 в точку (0.5, 0.1, 0.1, 0.1).
        Расстояние должно зависеть только от координаты X.
        """
        s = np.array([1.7, 0.1, 0.1, 0.1])
        
        # Расстояние (1.7 - 0.5) = 1.2
        expected_dist = 1.2
        
        result = dist_to_s(self.vor, s, self.vor.max_len)
        
        self.assertAlmostEqual(result, expected_dist, places=5, 
                               msg="Проекция на грань работает некорректно")
        





if __name__ == '__main__':
    unittest.main()