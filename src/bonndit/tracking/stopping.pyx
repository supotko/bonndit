#%%cython --annotate
#cython: language_level=3, boundscheck=False, wraparound=False, warn.unused=True, warn.unused_args=True,
# warn.unused_results=True
import os

from bonndit.utilc.cython_helpers cimport sub_vectors, angle_deg, sum_c, set_zero_matrix, bigger, smaller, mult_with_scalar, norm
import numpy as np
from bonndit.tracking.ItoW cimport Trafo
DTYPE = np.float64



cdef class Validator:
	def __cinit__(self, double[:,:,:] wm_mask, int[:] shape, double min_wm, str inclusion, double max_angle, Trafo trafo):
		self.min_wm = min_wm
		self.wm_mask = wm_mask
		self.shape = shape
		if inclusion:
			print(1)
			self.ROI = ROIValidator(inclusion)
		else:
			self.ROI = ROINotValidator(inclusion)
		if max_angle > 0:
			print(1)
			self.Curve = CurvatureValidator(max_angle, trafo)
		else:
			self.Curve = CurvatureNotValidator(max_angle, trafo)



	cdef bint wm_checker(self, double[:] point) nogil except *:
		""" Checks if the wm density is at a given point below a threshold.
		@param point: 3 dimensional point
		"""
		if self.wm_mask[int(point[0]), int(point[1]), int(point[2])] < self.min_wm:
			return True
		else:
			return False

	cdef bint index_checker(self, double[:] point) nogil except *:
		"""
		Checks if the index is within the array.
		@param point: 3 dimensional point
		@return: True if the point is not valid.
		"""
		if point[0] < 0 or point[1] < 0 or point[2] < 0:
			return True
		elif point[0] >= self.shape[0] or point[1] >= self.shape[1] or point[2] >= self.shape[2]:
			return True
		else:
			return False



	cdef bint next_point_checker(self, double[:] point) nogil except *:
		"""
		Check if a given direction is valid e.g. not zero and not infinity
		@param point: given direction
		@return:
		"""
		if sum_c(point) == 0 or sum_c(point) != sum_c(point):
			return True
		else:
			return False

	cdef void set_path_zero(self, double[:,:] path, double[:,:] features) nogil except *:
		set_zero_matrix(path)
		set_zero_matrix(features)


cdef class CurvatureNotValidator:
	def __cinit__(self, max_angle, trafo):
		self.max_angle = max_angle
		self.angle = 0
		self.points = np.zeros([5,3])
		self.trafo = trafo

	cdef bint curvature_checker(self, double[:,:] path,  double[:] features) nogil except *:
		return False

cdef class CurvatureValidator(CurvatureNotValidator):
	#def __cinit__(self, double max_angle):
#		super().__cinit__(max_angle)

	cdef bint curvature_checker(self, double[:,:] path,  double[:] features) nogil except *:
			"""
			Checks the angles between the current direction and the directions anlong the polygon. If a angle is to large returns True
			@param path: polygon to check
			@param features: save the angle between the current direction and the direction k points ago into the features.
			@return:
			"""
			cdef int l = 1, k = path.shape[0]
			cdef double length = 0
			self.trafo.itow(path[k])
			mult_with_scalar(self.points[3], 1, self.trafo.point_itow)
			sub_vectors(self.points[1], path[k-1], self.points[3])
			length += norm(self.points[1])
			while k >= 2 and length < 30 and l < k:
				l += 1
				sub_vectors(self.points[0], path[k-l], path[k-l+1])
				length += norm(self.points[0])
				self.angle = angle_deg(self.points[1], self.points[0])
				if self.angle > self.max_angle:
					return True
			else:
				features[0] = self.angle
				return False

cdef class ROINotValidator:
	def __cinit__(self, str inclusion):
		self.inclusion = np.zeros([3,3])
		self.inclusion_num = 0
		self.inclusion_check = np.zeros(1)


	cdef void included(self, double[:] point) nogil except *:
		pass

	cdef bint included_checker(self) nogil except *:
		return False

cdef class ROIValidator(ROINotValidator):
	def __cinit__(self, str inclusion):
		cubes = [os.path.join(inclusion, x) for x in os.listdir(inclusion) if x.endswith('.pts')]
		output = np.zeros((len(cubes)*2, 3))
		for i, cube in enumerate(cubes):
			points = open(cube)
			points = np.array([list(map(float, point.split())) for point in points])
			points = np.vstack((np.min(points, axis=0), np.max(points, axis=0)))
			output[2*i:2*(i+1)] = points
		self.inclusion = output
		self.inclusion_num = len(cubes)
		self.inclusion_check = np.zeros(len(cubes))

	cdef void included(self, double[:] point) nogil except *:
		cdef int i
		if sum_c(self.inclusion_check) == self.inclusion_num:
			return

		for i in range(self.inclusion_num):
			if smaller(point, self.inclusion[2*i]):
				continue
			if bigger(point, self.inclusion[2*i + 1]):
				continue
			self.inclusion_check[i] = 1



	cdef bint included_checker(self) nogil except *:
		return sum_c(self.inclusion_check) != self.inclusion_num
