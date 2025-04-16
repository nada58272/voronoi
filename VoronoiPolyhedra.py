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


class VoronoiPolyhedra(Voronoi):
	def __init__(self, grid):
		self.grid = grid
  
		self.coords4 = [] # координаты центров
		for var in itertools.product(range(-2, 2), repeat=4):
			self.coords4.append((self.grid.T).dot(var).tolist())
   
		# строим диаграмму Вороного
		super().__init__(self.coords4)
    
	def fillData(self):
		#! переименовать в central_point_index
		self.centr_index = self.coords4.index([0.0, 0.0, 0.0, 0.0])

	def findCentral(self):
		# находим суммарные расстояния от вершин многогранников до (0, 0, 0, 0) (sum_dist), затем находим минимальное 
		# значение. Таким образом находим ближайший центральный многогранник central

		sum_dist_min = 1000

		for region in self.regions:
			
			len_r = len(region)

			if -1 in region or len(region) == 0: continue # регион неограничен
				
			length = 0
			
			# для каждого региона перебираем все его вершины и для каждой ищем расстояние до [0.0, 0.0, 0.0, 0.0]
			for j in range(len_r):
				length += distance.euclidean(self.vertices[region[j]], [0.0, 0.0, 0.0, 0.0])

			if length < sum_dist_min:
				sum_dist_min = length
				v_min = region

		# ищу индекс минимального региона
		v_min = next((i for i, arr in enumerate(self.regions) if np.array_equal(arr, v_min)), None)

		print ('суммарное расстояние =', sum_dist_min, 'индекс центрального региона =', v_min)

		#! переименовать в central_region_index
		self.central_index = self.regions[v_min] # индексы координат (в coord4) центрального региона
		self.central = self.vertices[self.regions[v_min]]

	def findFaces3d(self):
		self.faces3d = []
		self.central_3d = []
		...
