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


from scipy.spatial import distance

from VoronoiPolyhedra import *

from FindSubGrids import *




# через минимальное расстояние до вершины

def dist_to_s(polyhedrons, s, vor4, delaunay, max_len):
    
    min_dist_to_pol = float('inf')
    index_proj = -1
    min_dist_vert_to_s = float('inf')
    min_vert = -1

    # находим ближайшую вершину к точке s
    
    for index in range(len(vor4.central)):

        dist_vert_to_s = distance.euclidean(s, vor4.central[index])

        if dist_vert_to_s < min_dist_vert_to_s:
            min_dist_vert_to_s = dist_vert_to_s
            min_vert = index
        

    # находим расстояние и проекцию на центральный многогранник
    
    for i in vor4.vertex_to_faces[min_vert]: # рассматриваем только грани,которым принадлежит ближайшая в точке s вершина
        # coord0 = s - (d0 + polyhedrons[i].bias)* polyhedrons[i].normal
        d0 = polyhedrons[i].normal @ (s - polyhedrons[i].center)
        coord0 = s - d0 * polyhedrons[i].normal
        simplex = delaunay.find_simplex(coord0)
        
        if simplex != -1: # если проекция принадлежит центральному многораннику
            dist = abs(d0)
            coords_to_central = coord0

            if min_dist_to_pol > dist:
                min_dist_to_pol = dist
                coords_proj = coords_to_central
                index_proj = i

        else:
            for face2d in polyhedrons[i].faces: #cycle

                d1 = face2d.normal @ (coord0 - face2d.center)
                # coord1 = coord0 - (d1 + face2d.bias)* polyhedrons[i].normal
                coord1 = coord0 - d1 * face2d.normal
                simplex = delaunay.find_simplex(coord1)

                if simplex != -1: # если проекция принадлежит центральному многораннику
                    dist = distance.euclidean(s, coord1)#, dtype = 'float')
                    coords_to_central = coord1

                    if min_dist_to_pol > dist:
                        min_dist_to_pol = dist
                        coords_proj = coords_to_central
                        index_proj = i

                else:

                    for edge in face2d.edges:
                        d2 = edge.normal @ (coord1 - edge.center)
                        #coord2 = coord1 - (d2 + edge.bias) * edge.normal
                        coord2 = coord1 - d2 * edge.normal
                        simplex = delaunay.find_simplex(coord2)

                        if simplex != -1: # если проекция принадлежит центральному многораннику
                            dist = distance.euclidean(s, coord2)#, dtype = 'float')
                            coords_to_central = coord2

                            if min_dist_to_pol > dist:
                                min_dist_to_pol = dist
                                coords_proj = coords_to_central
                                index_proj = i

                        else:
                            d3 = distance.euclidean(s, edge.vertex1)
                            d4 = distance.euclidean(s, edge.vertex2)

                            if d3 < d4:
                                dist = d3
                                coords_to_central = edge.vertex1
                            else:
                                dist = d4
                                coords_to_central = edge.vertex2

                            if min_dist_to_pol > dist:
            
                                min_dist_to_pol = dist
                                coords_proj = coords_to_central
                                index_proj = i
        
        
        # если расстояние до какой-либо грани < 1, то дальше не считаем
        if  min_dist_to_pol < 1:

            return min_dist_to_pol, coords_proj, index_proj

        
    return min_dist_to_pol * 2 / max_len, coords_proj, index_proj



# ищет точку s для подрешетки

def s_point(sub_grid, vor4):
    digits = np.array(list(product(range(-2, 2), repeat=4)))  # Все комбинации
    coords4_1 = np.dot(digits, sub_grid)  # Все точки подрешётки
    
    # Удаляем нулевую точку
    coords4_1 = coords4_1[~np.all(coords4_1 == 0, axis=1)]
    
    # Минимальное расстояние между точками
    if len(coords4_1) >= 2:
        dist_matrix = distance.pdist(coords4_1)
        min_dist = np.min(dist_matrix)
        if min_dist <= 2:
            return np.array([0, 0, 0, 0])


    # Поиск общей точки с coords4
    common_coord = np.array([x for x in coords4_1 if x.tolist() in vor4.coords4])
    if len(common_coord) == 0:
        return np.array([0, 0, 0, 0])
    

    # Ближайшая к центру точка
    dist_to_center = np.linalg.norm(common_coord, axis=1)
    closest_idx = np.argmin(dist_to_center)
    return common_coord[closest_idx] * 0.5