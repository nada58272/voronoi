from Polyhedra import *
from Factorization import *

#--------------------------------------------------------------------------------
def update_min_distance(dist, coords, index,
                        current_min_dist, current_coords, current_index):
    """
    Обновляет текущее минимальное расстояние и соответствующие координаты и индекс,
    если новое расстояние меньше.

    :param dist: float — новое расстояние
    :param coords: np.array — координаты точки с этим расстоянием
    :param index: int — индекс проекции (например, индекс многогранника)
    :param current_min_dist: float — текущее минимальное расстояние
    :param current_coords: np.array — текущие координаты проекции
    :param current_index: int — текущий индекс

    :return: (float, np.array, int) — обновлённые значения (min_dist, coords_proj, index_proj)
    """
    if dist < current_min_dist:
        return dist, coords, index
    else:
        return current_min_dist, current_coords, current_index

#--------------------------------------------------------------------------------

def find_faces_from_nearest_vertices(s, central, vertex_to_faces):
    """
    Находит грани, связанные с ближайшими вершинами к заданной точке s.

    :param s: np.array, координаты точки, от которой ищется расстояние.
    :param central: list of np.array, список вершин.
    :param vertex_to_faces: dict, словарь, где ключ - индекс вершины, значение - список граней, содержащих эту вершину.
    :return: set, множество общих граней, связанных с ближайшими вершинами.
    """
    min_dist_vert_to_s = float('inf') # минимальное расстояние до точки s
    min_dist_to_s_list = [] # список индексов ближайших вершин

    # Находим ближайшие вершины к точке s
    for index in range(len(central)):
        # считаем расстояние от вершины до точки s
        dist_vert_to_s = distance.euclidean(s, central[index])

        # если расстояние равно минимальному, добавляем индекс в список
        if dist_vert_to_s == min_dist_vert_to_s:
            min_dist_to_s_list.append(index)

        # если расстояние меньше минимального, обновляем минимальное и очищаем список
        elif dist_vert_to_s < min_dist_vert_to_s:
            min_dist_vert_to_s = dist_vert_to_s
            min_dist_to_s_list = [index]

    # Находим пересечение граней, связанных с ближайшими вершинами
    selected_lists = [set(vertex_to_faces[i]) for i in min_dist_to_s_list]
    common_elements = set.intersection(*selected_lists)

    return common_elements

#--------------------------------------------------------------------------------

def dist_to_s(polyhedrons, s, vor4, max_len):
    
    min_dist_to_pol = float('inf') #минимальное расстояние до центрального многогранника
    index_proj = -1
    min_dist_vert_to_s = float('inf')
    min_vert = -1
    coords_proj = np.array([1, 1, 1, 1])

    # находим ближайшую вершину к точке s
    '''
    for index in range(len(vor4.central)):

        dist_vert_to_s = distance.euclidean(s, vor4.central[index])

        if dist_vert_to_s < min_dist_vert_to_s:
            min_dist_vert_to_s = dist_vert_to_s
            min_vert = index
    ''' 
    # находим список ближайших вершин к точке s   
    nearest_faces = find_faces_from_nearest_vertices(s, vor4.central, vor4.vertex_to_faces)

    # находим расстояние и проекцию на центральный многогранник
    
    for i in nearest_faces:#vor4.vertex_to_faces[min_vert]: # рассматриваем только грани,которым принадлежит ближайшая в точке s вершина
        # coord0 = s - (d0 + polyhedrons[i].bias)* polyhedrons[i].normal
        d0 = polyhedrons[i].normal @ (s - polyhedrons[i].center)
        coord0 = s - d0 * polyhedrons[i].normal
        simplex = vor4.delaunay.find_simplex(coord0)
        
        if simplex != -1: # если проекция принадлежит центральному многораннику
            dist = abs(d0)
            coords_to_central = coord0
            
            min_dist_to_pol, coords_proj, index_proj = update_min_distance(dist, 
                                                                            coords_to_central, 
                                                                            i,
                                                                            min_dist_to_pol,
                                                                            coords_proj,
                                                                            index_proj
                                                                            )


            #coords_to_central = s - (d0 + polyhedrons[i].bias) * polyhedrons[i].normal
        

        else:
            for face2d in polyhedrons[i].faces: #cycle

                d1 = face2d.normal @ (coord0 - face2d.center)
                # coord1 = coord0 - (d1 + face2d.bias)* polyhedrons[i].normal
                coord1 = coord0 - d1 * face2d.normal
                simplex = vor4.delaunay.find_simplex(coord1)
                
                if simplex != -1: # если проекция принадлежит центральному многораннику
                    dist = distance.euclidean(s, coord1)#, dtype = 'float')
                    coords_to_central = coord1
                    
                    min_dist_to_pol, coords_proj, index_proj = update_min_distance(dist, 
                                                                           coords_to_central, 
                                                                           i,
                                                                           min_dist_to_pol,
                                                                           coords_proj,
                                                                           index_proj
                                                                          )
                    
                    
                else:

                    for edge in face2d.edges:
                        d2 = edge.normal @ (coord1 - edge.center)
                        #coord2 = coord1 - (d2 + edge.bias) * edge.normal
                        coord2 = coord1 - d2 * edge.normal
                        simplex = vor4.delaunay.find_simplex(coord2)
                        
                        if simplex != -1: # если проекция принадлежит центральному многораннику
                            dist = distance.euclidean(s, coord2)#, dtype = 'float')
                            coords_to_central = coord2

                            
                            min_dist_to_pol, coords_proj, index_proj = update_min_distance(dist, 
                                                                           coords_to_central, 
                                                                           i,
                                                                           min_dist_to_pol,
                                                                           coords_proj,
                                                                           index_proj
                                                                          )

                        else:
                            d3 = distance.euclidean(s, edge.vertex1)
                            d4 = distance.euclidean(s, edge.vertex2)
                            #print('3, 4', d3, d4)
                            if d3 < d4:
                                dist = d3
                                coords_to_central = edge.vertex1
                            else:
                                dist = d4
                                coords_to_central = edge.vertex2

                            min_dist_to_pol, coords_proj, index_proj = update_min_distance(
                                                                        dist, 
                                                                        coords_to_central, 
                                                                        i,
                                                                        min_dist_to_pol,
                                                                        coords_proj,
                                                                        index_proj
                                                                        )
     
                            
        # если расстояние до какой-либо грани < 1, то дальше не считаем
    
        min_dist_to_pol = min_dist_to_pol * 2 / max_len
        
        if  min_dist_to_pol < 1:
            return min_dist_to_pol, coords_proj, index_proj

    return min_dist_to_pol, coords_proj, index_proj

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
