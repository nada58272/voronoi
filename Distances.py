from Polyhedra import *
from Factorization import *
import math

CHECK_DIST = True  # проверка расстояний по теореме Пифагора

#--------------------------------------------------------------------------------

def find_faces_from_nearest_vertices(s, central, vertex_to_faces):
    """
    Находит грани, связанные с ближайшими вершинами к заданной точке s.

    :param s: np.array, координаты точки, от которой ищется расстояние.
    :param central: list of np.array, список вершин.
    :param vertex_to_faces: массив, где индекс - индекс вершины в central, значение - список граней, содержащих эту вершину.
    :return: set, множество общих граней, связанных с ближайшими вершинами.
    """
    min_dist_vert_to_s = float('inf') # минимальное расстояние до точки s
    min_dist_to_s_list = [] # список индексов ближайших вершин

    # Находим ближайшие вершины к точке s
    for index in range(len(central)):
        # считаем расстояние от вершины до точки s
        dist_vert_to_s = distance.euclidean(s, central[index])

        if np.allclose((dist_vert_to_s - min_dist_vert_to_s), 0, atol=1e-7):
            min_dist_to_s_list.append(index)

        # если расстояние меньше минимального, обновляем минимальное и очищаем список
        else:
            
            if dist_vert_to_s < min_dist_vert_to_s:
                min_dist_vert_to_s = dist_vert_to_s
                min_dist_to_s_list = [index]

        # если расстояние равно минимальному, добавляем индекс в список
       # if dist_vert_to_s == min_dist_vert_to_s:
        #    min_dist_to_s_list.append(index)

        # если расстояние меньше минимального, обновляем минимальное и очищаем список
        #elif dist_vert_to_s < min_dist_vert_to_s:
         #   min_dist_vert_to_s = dist_vert_to_s
          #  min_dist_to_s_list = [index]

    # Находим пересечение граней, связанных с ближайшими вершинами
    selected_lists = [set(vertex_to_faces[i]) for i in min_dist_to_s_list]
    common_elements = set.intersection(*selected_lists)

    return common_elements

#--------------------------------------------------------------------------------

def checkDist(dist1, dist2):
    if not math.isclose(dist1, dist2, abs_tol=1e-9):
        print('dist_test', math.isclose(dist1, dist2, abs_tol=1e-9), dist1 - dist2)
        
#--------------------------------------------------------------------------------
    
def dist_to_s(vor4, s, max_len):
    polyhedrons = vor4.polyhedrons
    
    min_dist_to_pol = float('inf') #минимальное расстояние до центрального многогранника
    index_proj = -1
    coords_proj = np.array([1, 1, 1, 1])
    list_parts = []
 

    # находим список ближайших вершин к точке s   
    nearest_faces = find_faces_from_nearest_vertices(s, vor4.central, vor4.vertex_to_faces)
    
    def update_min_distance():
        nonlocal dist, coords_to_central, i#, list_parts_0
        nonlocal min_dist_to_pol, coords_proj, index_proj#, list_parts
        
        if dist < min_dist_to_pol:
            min_dist_to_pol = dist
            coords_proj = coords_to_central
            index_proj = i
            #list_parts = list_parts_0.copy() # сохраняем список граней, к которым принадлежит ближайшая вершина
 
    # находим расстояние и проекцию на центральный многогранник
    
    for i in nearest_faces: # рассматриваем только грани,которым принадлежит ближайшая в точке s вершина
        d0 = polyhedrons[i].normal @ (s - polyhedrons[i].center)
        coord0 = s - d0 * polyhedrons[i].normal
        simplex = vor4.delaunay.find_simplex(coord0, tol=1e-9) # находим проекцию на центральный многогранник
     
        d0_squared = d0 * d0
         
        
        if simplex != -1: # если проекция принадлежит центральному многораннику
            dist = abs(d0)
            coords_to_central = coord0
            #list_parts_0 = [i] # сохраняем список граней, к которым принадлежит ближайшая вершина
            #print(99)
            update_min_distance()
            continue     

        for face2d in polyhedrons[i].faces: #cycle

            d1 = face2d.normal @ (coord0 - face2d.center)
            # coord1 = coord0 - (d1 + face2d.bias)* polyhedrons[i].normal
            coord1 = coord0 - d1 * face2d.normal
            simplex = vor4.delaunay.find_simplex(coord1, tol=1e-9) # находим проекцию на центральный многогранник

            d1_squared = d1 * d1 # храним хквадрат расстояния до 2х мерной грани
            
            if simplex != -1: # если проекция принадлежит центральному многораннику
                dist = distance.euclidean(s, coord1)#, dtype = 'float')
                coords_to_central = coord1

                #list_parts_0 = [i, face2d]


                # проверяем, что посчитанное расстояние совпадает с расстоянием по теореме Пифагора
                dist_test = dist * dist - d1_squared - d0_squared
                if not math.isclose(dist_test, 0., abs_tol=1e-9):
                    print('dist_test, d1', math.isclose(dist_test, 0., abs_tol=1e-9), dist_test)
                
                update_min_distance()
                continue


            for edge in face2d.edges:
                d2 = edge.normal @ (coord1 - edge.center)
                #coord2 = coord1 - (d2 + edge.bias) * edge.normal
                coord2 = coord1 - d2 * edge.normal
                simplex = vor4.delaunay.find_simplex(coord2, tol=1e-9) # находим проекцию на центральный многогранник

                d2_squared = d2 * d2 # храним квадрат расстояния до ребра
                
                if simplex != -1: # если проекция принадлежит центральному многораннику
                    dist = distance.euclidean(s, coord2)#, dtype = 'float')
                    coords_to_central = coord2

                    #list_parts_0 = [i, face2d, edge]

                    # проверяем, что посчитанное расстояние совпадает с расстоянием по теореме Пифагора
                    dist_test = dist * dist - d1_squared - d0_squared - d2_squared
                    if not math.isclose(dist_test, 0., abs_tol=1e-9):
                        print('dist_test, d2', math.isclose(dist_test, 0., abs_tol=1e-9), dist_test)

                    update_min_distance()

                else:
                    d3 = distance.euclidean(coord2, edge.vertex1)                    
                    d4 = distance.euclidean(coord2, edge.vertex2)

                    if d3 < d4:
                        dist = distance.euclidean(s, edge.vertex1)
                        coords_to_central = edge.vertex1
                        d34_squared = d3 * d3 # храним квадрат расстояния до вершины 1
                        #list_parts_0 = [i, face2d, edge, edge.vertex1]

                    else:
                        dist = distance.euclidean(s, edge.vertex2)
                        coords_to_central = edge.vertex2
                        d34_squared = d4 * d4 # храним квадрат расстояния до вершины 1
                        #list_parts_0 = [i, face2d, edge, edge.vertex2]

                    # проверяем, что посчитанное расстояние совпадает с расстоянием по теореме Пифагора
                    if CHECK_DIST:
                        checkDist(dist * dist, d0_squared + d1_squared + d2_squared + d34_squared)

                    update_min_distance()
    
                            
        # если расстояние до какой-либо грани < 1, то дальше не считаем
    
        #min_dist_to_pol = min_dist_to_pol * 2 / max_len
        
        if  min_dist_to_pol * 2 / max_len < 1:
            return min_dist_to_pol * 2 / max_len
        
        
        
    min_dist_to_pol = min_dist_to_pol * 2 / max_len

    return min_dist_to_pol#, coords_proj, index_proj, list_parts

#--------------------------------------------------------------------------------
# Достаточный набор комбинаций

digits = np.array(list(product(
	[0, 1, 2], 
	[0, 1, -1, 2, -2], 
	[0, 1, -1, 2, -2], 
	[0, 1, -1, 2, -2], 
)))

# Удаляем нулевую точку
#! можно оптимизировать еще на этом этапе: не включать в список элементы вида [0, negative, ..., ...], [0, 0, negative, ...], [0, 0, 0, negative]
digits = digits[~np.all(digits == 0, axis=1)]

#--------------------------------------------------------------------------------
# ищет точку s для подрешетки
# digits уже определен!

def s_point(sub_grid, vor4):   
    coords4_1 = np.dot(digits, sub_grid)  # Все точки подрешётки
    
    # Поиск общей точки с coords4
    common_coord = np.array([x for x in coords4_1 if x.tolist() in vor4.coords4])
    if len(common_coord) == 0:
        return np.array([0, 0, 0, 0])
    
    # Ближайшая к центру точка
    dist_to_center = np.linalg.norm(common_coord, axis=1)
    closest_idx = np.argmin(dist_to_center)
    return common_coord[closest_idx] * 0.5

#--------------------------------------------------------------------------------

def center_points(sub_grid):       
    grid_points = np.dot(digits, sub_grid)  # Все точки подрешётки для данного набора комбинаций
    points = []
    
    # оцениваем минимальное расстояние между центрами многогранников
    central_dist_min = np.linalg.norm(grid_points[0])
    for point in grid_points:
        central_dist = np.linalg.norm(point)
        if central_dist < central_dist_min:
            central_dist = central_dist_min

    # отбираем точки для тех многогранников, расстояние от центра до центра которых не более central_dist_min + diameter
    for point in grid_points:
        central_dist = np.linalg.norm(point)
        if central_dist < central_dist_min + 2:
            points.append(np.array(point))

    return points
