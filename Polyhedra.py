import numpy as np
from scipy.spatial import Voronoi, Delaunay, ConvexHull, distance
from itertools import product, combinations

#--------------------------------------------------------------------------------

class Edge2D:
    def __init__(self, vertex1, vertex2, face_center):
        '''
        Инициализация ребра.
        :vertex1: Координаты первой вершины ребра.
        :vertex2: Координаты второй вершины ребра.
        :face_center: Координаты центра грани.
        '''
        self.vertex1 = np.array(vertex1)
        self.vertex2 = np.array(vertex2)
        self.face_center = np.array(face_center)

        # Вычисляем середину ребра
        self.center = (self.vertex1 + self.vertex2) / 2

        # Вычисляем нормаль к ребру
        #edge_vector = self.vertex2 - self.vertex1
        self.normal = self.center - self.face_center  # Перпендикулярный вектор
        self.normal = self.normal / np.linalg.norm(self.normal)  # Нормализация
        
        self.bias = self.normal @ vertex1

        # Проверяем направление нормали (она должна быть направлена от центра грани)
        vector_to_center = self.face_center - self.center
        if np.dot(self.normal, vector_to_center) > 0:
            self.normal = -self.normal  # Меняем направление нормали

    def __repr__(self):
        return f"Edge2D(vertex1={self.vertex1}, vertex2={self.vertex2}, normal={self.normal}, \
                midpoint={self.center})"

#--------------------------------------------------------------------------------

class Face2D:
    def __init__(self, vertices, polyhedron_center, vor4):
        '''
        Инициализация 2D грани многогранника.
        '''
        self.vertices = vertices  # Координаты вершин грани
        self.vor4 = vor4
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
        self.normal = np.array(self.center) - np.array(self.parent_center)

        # Нормализуем нормаль
        self.normal = self.normal / np.linalg.norm(self.normal)
    
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
            self.edges.append(Edge2D(vertex1, vertex2, self.center))

    
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
            face = Face2D(face_vertices, self.center, self.vor4)
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
				#print(common_coords, y[comb[0]], y[comb[1]])
				
		

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

		self.list_neigh_points = [] # список центров многогранников, граничащих с центральным
		self.list_neigh_points_ind = [] # индексы из coords4

		for pair in self.ridge_points:

			if self.central_point_index not in pair: continue

			if pair[0] == self.central_point_index:
				
				# если в списке 4 многогранника есть вершины хоть одой 3 грани центрального многогранника
				# то этот многогранник соседний с центральным
				if any(all(item in self.regions[self.point_region[pair[1]]] for item in small) for small in self.faces_3d):
					self.list_neigh_points.append(self.coords4[pair[1]])
					self.list_neigh_points_ind.append(pair[1])
				
			else:
				if any(all(item in self.regions[self.point_region[pair[0]]] for item in small) for small in self.faces_3d):
					self.list_neigh_points.append(self.coords4[pair[0]])
					self.list_neigh_points_ind.append(pair[0])

		# ищем для каждого центра ближайшую трехмерную грань
		# берем каждую точку из списка соседних центральных точек и ищем расстояние от нее до каждой 3х мерной грани
		# суммируем и берем минимальное ссумарное расстояние. это и будет 3х мерная грань между точкой (0, 0, 0, 0)
		# точкой из списка
		# и вектор, соединяющий центры будет вектором нормали к данной 3х мерной грани

		self.list_ridge_edge = [] # список граней, в том же порядке, что и соседние точки (list_neigh_points)


		# перебираем все соседние центры
		for point in self.list_neigh_points:
			
			min_dist_point_edge = 1000

			# для каждой 3х мерной грани   
			for edge in self.edge_central_coords:
				
				dist_point_edge = 0
				
				# для каждой вершины их 3х мерной грани считаем растояние до соседнего центра и суммируем
				for coord in edge:
					
					dist_point_edge += distance.euclidean(coord, point)
					
				# если суммарное растояние меньше минимального, записываем новое суммарное растояние и индекс грани    
				if min_dist_point_edge > dist_point_edge:
					min_dist_point_edge = dist_point_edge

			self.list_ridge_edge.append(edge)


	# строим список 2d граней 3d граней
	def find_2d_subfaces(self):
		self.list_faces = []
		for face3d in self.faces_2d:
      
			list_vert_3d = []
			for face2d in face3d:
       
				list_verts = []
				for vert in face2d:
					list_verts.append(self.vertices[vert])

				list_vert_3d.append(list_verts) 

			self.list_faces.append(list_vert_3d)
		


	# нормируем вектора нормали
	# ищем точки, принадлежащие граням и находящиеся на векторах нормали
	def count_v_norm(self):
		self.v_norm = [] # здесь будем хранить нормированные вектора нормали
		self.point_face_norm = [] #точки, принадлежащие граням и находящиеся на векторах нормали

		for vec in self.list_neigh_points:
			self.v_norm.append(np.array(vec / np.linalg.norm(vec)))
			self.point_face_norm = np.array(vec) / 2

	def get_vertex_faces_3d(self):
		self.vertex_to_faces = [] # собой список списков 3D-граней для каждой вершины 
								  # (где индексы граней берутся из edge_central_coords)
		'''    
		берем каждую верщину из central и перебираем все 3х мерные грани - смотрим есть ли данная вершина в этой 
		3х мерной грани если да, то тдобавляем ее в список
		'''
		for vertex in self.central:
			
			tamp_list = [] # список граней для текущей вершины
			for face in range(len(self.edge_central_coords)):
				if np.any(np.all(self.edge_central_coords[face] == vertex, axis=1)):
					tamp_list.append(face)
					
			self.vertex_to_faces.append(tamp_list)

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
		self.findCentral()
		self.fillData()
		self.findFaces3d()
		self.findFaces2d()
		self.findEdges()
		self.findNeigh()

		# диаметр многогранника (ищем максимальное расстояние между вершинами) - 2
		self.max_len = round(2 * distance.euclidean(np.array([0, 0, 0, 0]), self.central[1]), 10)
		self.min_d = 1

		self.find_2d_subfaces()
		self.count_v_norm()
		self.createPolyhedrons()
		self.get_vertex_faces_3d()
  
		self.createTriangulation()
