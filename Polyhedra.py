import numpy as np
from scipy.spatial import Voronoi, Delaunay, ConvexHull, distance
from itertools import product, combinations

#--------------------------------------------------------------------------------

class Edge2D:
    def __init__(self, vertex1, vertex2, face_center, normal_2d, normal_3d):
        '''
        Инициализация ребра.
        :vertex1: Координаты первой вершины ребра.
        :vertex2: Координаты второй вершины ребра.
        :face_center: Координаты центра грани.
        '''

        self.vertex1 = np.array(vertex1) # координаты вершин
        self.vertex2 = np.array(vertex2) # координаты вершин
        self.face_center = np.array(face_center) # Центр 2D-грани
        self.normal_2d = normal_2d # Нормаль к 2D-грани
        self.normal_3d = normal_3d # Нормаль к 3D-гиперграни

        # Вычисляем середину ребра
        self.center = (self.vertex1 + self.vertex2) / 2

        # Вычисляем нормаль к ребру
        #edge_vector = self.vertex2 - self.vertex1
        #self.normal = self.center - self.face_center  # Перпендикулярный вектор
        #self.normal = self.normal / np.linalg.norm(self.normal)  # Нормализация
        normal_candidate = self._find_edge_normal(self.vertex1, self.vertex2, self.normal_2d, self.normal_3d)
        
        # Проверка на вырожденность и нормализация
        norm_val = np.linalg.norm(normal_candidate)

        if norm_val < 1e-9:
            # print(f"Warning: Degenerate normal calculated for Edge2D ({self.vertex1}, {self.vertex2}). Setting to zero vector.")
            self.normal = np.zeros_like(normal_candidate)
            self.bias = 0.0
            
        else:
            oriented_normal = normal_candidate / norm_val

            # Проверяем и корректируем направление нормали (она должна быть направлена "наружу" от центра 2D-грани)
            vector_to_face_center = self.face_center - self.center
            if np.dot(oriented_normal, vector_to_face_center) > 1e-7: # Допуск для численной стабильности
                self.normal = -oriented_normal
            else:
                self.normal = oriented_normal
            
            self.bias = self.normal @ self.vertex1

        #self.bias = self.normal @ self.vertex1

        # Проверяем направление нормали (она должна быть направлена от центра грани)
        #vector_to_center = self.face_center - self.center
        #if np.dot(self.normal, vector_to_center) > 0:
         #   self.normal = -self.normal  # Меняем направление нормали
    
    def _find_edge_normal(self, vertex1, vertex2, n_2d, n_3d):
            # Вектор ребра
            edge = vertex2 - vertex1
            
            # Проверка, не вырождено ли само ребро
            if np.linalg.norm(edge) < 1e-9:
                # print(f"Warning: Degenerate edge in Edge2D ({v1}, {v2}). Returning zero normal.")
                return np.zeros_like(vertex1) # Возвращаем нулевой вектор той же размерности
            
            edge = edge / np.linalg.norm(edge)
   # Матрица условий ортогональности
            #A = np.vstack([edge, n_2d, n_3d]).T
            A = np.column_stack((edge, n_2d, n_3d))

            # SVD разложение

            try:
                U, S, Vt = np.linalg.svd(A)
                
                if S[-1] < 1e-7: # Порог для определения вырожденности/линейной зависимости
                  return np.zeros_like(vertex1) # 
    
    # Нормаль к ребру - последний столбец U
                n_edge_candidate = U[:, -1]

            except np.linalg.LinAlgError:
                # print(f"Warning: SVD failed in Edge2D for edge ({v1}, {v2}). Returning zero normal.")
                return np.zeros_like(vertex1)
   
            # Нормировка
            #n_edge = n_edge / np.linalg.norm(n_edge)

            return n_edge_candidate

    def __repr__(self):
        return f"Edge2D(vertex1={self.vertex1}, vertex2={self.vertex2}, normal={self.normal}, \
                midpoint={self.center})"


#--------------------------------------------------------------------------------

class Face2D:
    def __init__(self, vertices, polyhedron_center, vor4, normal_parent):
        '''
        Инициализация 2D грани многогранника.
        '''
        self.vertices = vertices  # Координаты вершин грани
        self.vor4 = vor4
        self.parent_normal = normal_parent # нормаль к многограннику
        self.parent_center = polyhedron_center # центр 3х мерной грани
        self._calculate_center()  # Центр грани
        self._calculate_normal()  # Нормаль к грани
        self._calculate_edges()  # ребра в 2х мерном многограннике
        self.bias = self.normal @ vertices[0]

    def _calculate_center(self):
        '''
        Вычисляет центр грани как среднее арифметическое координат вершин.
        '''
        self.center = tuple(np.mean(self.vertices, axis=0))

    def _calculate_normal(self):
        # Вычисляем нормаль к грани как вектор из центра 3х мерного многогранника в центр 2х мерной грани
        #self.normal = np.array(self.center) - np.array(self.parent_center)
        v1 = self.vertices[0] - self.vertices[1]
        v2 = self.vertices[0] - self.vertices[2]
        A = np.column_stack((v1, v2, self.parent_normal))
        U, S, Vt = np.linalg.svd(A)
        # который ортогонален пространству столбцов A.
        normal_candidate = U[:, -1] # Используем более нейтральное имя
		
		# Нормализация
        norm_val = np.linalg.norm(normal_candidate)
        if norm_val < 1e-9: # Защита от деления на ноль
			# print(f"Warning: Degenerate normal calculated in Face2D for vertices {self.vertices}. Setting to zero vector.")
            self.normal = np.zeros_like(normal_candidate)
            self.bias = 0.0
            return
			
        oriented_normal = normal_candidate / norm_val # Имя для нормализованной, но еще не ориентированной
		
		# Проверка и коррекция направления
        vector_to_parent_center = np.array(self.parent_center) - np.array(self.center)
        if np.dot(oriented_normal, vector_to_parent_center) > 0:
			# Если скалярное произведение > 0, значит нормаль и вектор к центру родителя сонаправлены.
			# Это означает, что нормаль указывает "внутрь" относительно parent_center.
			# Меняем направление, чтобы она указывала "наружу".
            final_normal = -oriented_normal
        else:
            final_normal = oriented_normal
		
		# Обновляем атрибуты экземпляра!
        self.normal = final_normal
    
    # поиск ребер для 2х мерной грани (используется при создании класса)
    def find_edges(self):
        list_ind = [] # список вершин в индексах
        list_edges_face_2d_ind = [] # список ребер грани в индексах
        #list_edges_face_2d_coord = [] # список ребер грани в координатах
        
        for vert in self.vertices:
            ind = np.argwhere((self.vor4.vertices == vert).all(axis=1))[0][0]
            list_ind.append(ind)
            
        for edge in self.vor4.list_edges:
            if edge[0] in list_ind and edge[1] in list_ind:
                list_edges_face_2d_ind.append(edge)
                
        
        return self.vor4.vertices[list_edges_face_2d_ind]

    def _calculate_edges(self):
        # составляем список ребер (расстояние между вершинами, соединенными ребрами = 1)
        self.edges = []
        
        found_edges = self.find_edges()        
        
        for vertex1, vertex2 in found_edges:
			# поменяла 12.09.25 self.parent_center на self.parent_normal
			
            self.edges.append(Edge2D(vertex1, vertex2, self.center, self.normal, self.parent_normal))
            #self.edges.append(Edge2D(vertex1, vertex2, self.center))
		
    
    def __repr__(self):
        '''
        Возвращает строковое представление грани.
        '''
        return (f"Face2D(vertices={self.vertices}), center={self.center}, "
                f"normal={self.normal}")

#--------------------------------------------------------------------------------

class Face3D:
    def __init__(self, vertices, faces_list, normal, vor4):
        '''
        Инициализация многогранника в 4-мерном пространстве.
        :vertices: Список координат вершин многогранника .
        :normal: Вектор нормали к многограннику.
        '''

        self.vor4 = vor4
        self.vertices = vertices  # Координаты вершин
        self.normal = normal  # Нормаль к многограннику
        self.faces_list = faces_list  # 2d Грани 3d многогранника
        self._calculate_center()  # Центр многогранника
        self._find_faces() # грани многогранника (как экземпляр класса)
        self.bias = normal @ vertices[0] # смещение
    

    def _calculate_center(self):
        '''
        Вычисляет центр многогранника как среднее арифметическое координат вершин.
        :return: Координаты центра.
        '''
        self.center = tuple(np.mean(self.vertices, axis=0))

    
    def _find_faces(self):
        '''
        Находит грани многогранника и создаёт объекты Face2D.
        :return: Список объектов Face2D.
        '''
        
        self.faces = []

        for face_vertices in self.faces_list:
            face = Face2D(face_vertices, self.center, self.vor4, self.normal)
            self.faces.append(face)
            
        #for index in range(len(self.faces_list)):
            #face = Face2D(self.faces_list, index, self.center)
            #self.faces.append(face)
    

    def __repr__(self):
        '''
        Возвращает строковое представление многогранника.
        '''
        return (f"Face3D(vertices={self.vertices}, normal={self.normal}, "
                f"center={self.center}, faces={self.faces})")

#--------------------------------------------------------------------------------

class VoronoiPolyhedra(Voronoi):
	def __init__(self, grid):
		self.grid = grid
  
		self.coords4 = [] # координаты центров

		# меньше 3 нельзя - не хватает точке для подсчета
		for var in product(range(-3, 3), repeat=4):
			self.coords4.append((self.grid.T).dot(var).tolist())
   
		# строим диаграмму Вороного
		super().__init__(self.coords4)
    
	def fillData(self):
		self.central_point_index = self.coords4.index([0.0, 0.0, 0.0, 0.0])

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
				#print(sum_dist_min, v_min)

		# ищу индекс минимального региона
		v_min = next((i for i, arr in enumerate(self.regions) if np.array_equal(arr, v_min)), None)

		print ('суммарное расстояние =', sum_dist_min, 'индекс центрального региона =', v_min)

		self.central_region_index = self.regions[v_min] # индексы координат (в coord4) центрального региона
		self.central = self.vertices[self.regions[v_min]]

	def findFaces3d(self):
		
		self.faces_3d = [] # список 3d граней центрального региона
		self.central_3d = []
		self.to_remove = [] # список не 3d граней (список граней к удалению из зх мерных граней)

		# находим 3d грани

		for ridge in self.ridge_vertices:
			
			if np.all(np.isin(ridge, self.central_region_index)):
				self.central_3d.append(self.vertices[ridge])
				self.faces_3d.append(ridge)

		# в данный список могли попасть не только 3х мерные грани
		# если какая то грань полностью входит в другую, то это не 3х мерная грань

		for face_1 in self.faces_3d:

			for face_2 in self.faces_3d:
				# если face_1 входит в face_2 и длина face_1 меньше  face_2, то face_1 не 3х мерная грань
				if np.all(np.isin(face_1, face_2)) and len(face_1) < len(face_2):
					self.to_remove.append(face_1)
					break
						
		#  удаляем лишние грани
		for face in self.to_remove:    
			self.faces_3d.remove(face)


		# переводим в координаты индексы 3d грани

		self.edge_central_coords = []  # Координаты вершин 3D-граней

		for face in self.faces_3d:
			list_c = []
			for vert in face:
				list_c.append(self.vertices[vert])
			
			self.edge_central_coords.append(list_c)
		


	def findFaces2d(self):
	# формируем список 2d граней (если у 2х 3d граней есть общие вршины и из более 2, то это 2d грань)

		self.faces_2d  = [[] for _ in range(len(self.faces_3d))] # список 2d граней с разбивкой по 3d граням

		# берем все комбиации 3d граней по 2 и ищем общие вершины
		for comb in combinations(range(len(self.faces_3d)), 2):
			common_coords = list(set(self.faces_3d[comb[0]]) & set(self.faces_3d[comb[1]]))
			
			if len(common_coords) > 2:
				self.faces_2d[comb[0]].append(common_coords)
				self.faces_2d[comb[1]].append(common_coords)
				
		

	def findEdges(self):  
		self.edges  = [[] for _ in range(len(self.faces_2d))]

		for index in range(len(self.faces_2d)):
			
			for comb in combinations(range(len(self.faces_2d[index])), 2):
				common_coords = list(set(self.faces_2d[index][comb[0]]) & set(self.faces_2d[index][comb[1]]))
				
				if len(common_coords) == 2:
					common_coords.sort()
					self.edges[index].append(common_coords)

		self.list_edges = []

		for face in self.edges:
			for edge in face:
				if edge not in self.list_edges:
					self.list_edges.append(edge)


	def findNeigh(self):
		# составляем список центров многогранников, граничащих с центральным
		# индексы координат в coords4

		self.list_neigh_points = [[] for _ in range(len(self.faces_3d))] # список центров многогранников, граничащих с центральным
		self.list_neigh_points_ind = [[] for _ in range(len(self.faces_3d))] # индексы из coords4



		for pair in self.ridge_points:

			if self.central_point_index not in pair: continue

			if pair[0] == self.central_point_index:
				
				region = self.regions[self.point_region[pair[1]]]

				# если в списке 4 многогранника есть вершины хоть одой 3 грани центрального многогранника
				# то этот многогранник соседний с центральным
				for i in range(len(self.faces_3d)):
					if all(item in region for item in self.faces_3d[i]):
						self.list_neigh_points_ind[i] = pair[1]
						self.list_neigh_points[i] = self.coords4[pair[1]]


				#if any(all(item in self.regions[self.point_region[pair[1]]] for item in small) for small in self.faces_3d):
				#	self.list_neigh_points.append(self.coords4[pair[1]])
				#	self.list_neigh_points_ind.append(pair[1])
				
			else:
				region = self.regions[self.point_region[pair[0]]]
				for i in range(len(self.faces_3d)):
					if all(item in region for item in self.faces_3d[i]):
						self.list_neigh_points_ind[i] = pair[0]
						self.list_neigh_points[i] = self.coords4[pair[0]]

				#if any(all(item in self.regions[self.point_region[pair[0]]] for item in small) for small in self.faces_3d):
				#	self.list_neigh_points.append(self.coords4[pair[0]])
				#	self.list_neigh_points_ind.append(pair[0])

		# ищем для каждого центра ближайшую трехмерную грань
		# берем каждую точку из списка соседних центральных точек и ищем расстояние от нее до каждой 3х мерной грани
		# суммируем и берем минимальное ссумарное расстояние. это и будет 3х мерная грань между точкой (0, 0, 0, 0)
		# точкой из списка
		# и вектор, соединяющий центры будет вектором нормали к данной 3х мерной грани

		#self.list_ridge_edge = [] # список граней, в том же порядке, что и соседние точки (list_neigh_points)


		# перебираем все соседние центры
		#for point in self.list_neigh_points:
			
		#	min_dist_point_edge = 1000

			# для каждой 3х мерной грани   
		#	for edge in self.edge_central_coords:
				
		#		dist_point_edge = 0
				
				# для каждой вершины их 3х мерной грани считаем растояние до соседнего центра и суммируем
		#		for coord in edge:
					
		#			dist_point_edge += distance.euclidean(coord, point)
					
				# если суммарное растояние меньше минимального, записываем новое суммарное растояние и индекс грани    
		#		if min_dist_point_edge > dist_point_edge:
		#			min_dist_point_edge = dist_point_edge

		#	self.list_ridge_edge.append(edge)

    # строим список 2d граней 3d граней
	def find_2d_subfaces(self):
		self.list_faces = [] # список 2d граней 3d граней в координатах

		# перебираем все 3х мерные грани 
		# self.faces_2d - это список списков 2х мерных граней, которые принадлежат данной 3х мерной грани в индексах
		for face3d in self.faces_2d:
      
			list_vert_3d = [] # список 2d граней в координатах
			# перебираем все 2х мерные грани, которые принадлежат данной 3х мерной грани
			for face2d in face3d:
       
				list_verts = [] # список вершин 2х мерной грани в координатах
				# перебираем все вершины 2х мерной грани и записываем их координаты
				for vert in face2d:
					list_verts.append(self.vertices[vert])

				list_vert_3d.append(list_verts) 

			self.list_faces.append(list_vert_3d)
		


	# нормируем вектора нормали
	# ищем точки, принадлежащие граням и находящиеся на векторах нормали
	def count_v_norm(self):
		self.v_norm = [] # здесь будем хранить нормированные вектора нормали
		#self.point_face_norm = [] #точки, принадлежащие граням и находящиеся на векторах нормали

		for vec in self.list_neigh_points:
			self.v_norm.append(np.array(vec / np.linalg.norm(vec)))
			#self.point_face_norm = np.array(vec) / 2

	def get_vertex_faces_3d(self):
		self.vertex_to_faces = [] # собой список списков 3D-граней для каждой вершины 
								  # (где индексы граней берутся из edge_central_coords)
		'''    
		берем каждую верщину из central и перебираем все 3х мерные грани - смотрим есть ли данная вершина в этой 
		3х мерной грани если да, то тдобавляем ее в список
		'''
		#for vertex in self.central:
		for vertex in self.central_region_index:
			
			temp_list = [] # список граней для текущей вершины
			
			#for face in range(len(self.edge_central_coords)):
			for face in range(len(self.faces_3d)):
				
				# если вершина входит в грань, то добавляем ее в список
				#if np.any(np.all(np.isclose(self.edge_central_coords[face], vertex, atol=1e-8), axis=1)):
				if vertex in self.faces_3d[face]:
					temp_list.append(face)

			# добавляем список граней для текущей вершины в общий список	
			self.vertex_to_faces.append(temp_list)

	def createPolyhedrons(self):
		# расчитываем все грани центрального многогранника

		self.polyhedrons = []  # Список для хранения всех объектов

		for face_index in range(len(self.edge_central_coords)):
			polyhedron = Face3D(self.edge_central_coords[face_index], self.list_faces[face_index], self.v_norm[face_index], self)
			self.polyhedrons.append(polyhedron)

	def createTriangulation(self):
        # строим триангуляцию
		self.delaunay = Delaunay(self.central)

		# Получаем все вершины симплексов
		simplices_vertices = self.delaunay.simplices

		# Уникальные вершины
		unique_vertices = np.unique(simplices_vertices)

		# Проверяем, что все точки из central используются
		if len(unique_vertices) == len(self.central):
			print('Все точки используются в триангуляции.')
		else:
			print('Не все точки используются в триангуляции.')

		# Выпуклая оболочка исходных точек
		original_hull = ConvexHull(self.central)

		# Все точки триангуляции
		triangulated_points = self.central[self.delaunay.simplices.flatten()]

		# Выпуклая оболочка точек триангуляции
		triangulated_hull = ConvexHull(triangulated_points)

		# Сравнение объемов выпуклых оболочек
		is_convex_hull_correct = np.isclose(original_hull.volume, triangulated_hull.volume)
		print('Выпуклость триангуляции сохранена:', is_convex_hull_correct)

	def build(self):
		self.fillData()
		self.findCentral()
		#self.fillData()
		self.findFaces3d()
		self.findFaces2d()
		self.findEdges()
		self.findNeigh()

		def find_max_len():
			max_len = 0
			for i in range(len(self.central)):
				a = distance.euclidean(self.central[i], np.array([0, 0, 0, 0]))
				if a > max_len:
					max_len = a
			return max_len * 2
		
		# диаметр многогранника (ищем максимальное расстояние между вершинами) - 2
		self.max_len = find_max_len()
		self.min_d = 1

		self.find_2d_subfaces()
		self.count_v_norm()
		self.createPolyhedrons()
		self.get_vertex_faces_3d()
  
		self.createTriangulation()
