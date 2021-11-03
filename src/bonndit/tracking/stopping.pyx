#%%cython --annotate
#cython: language_level=3, boundscheck=True, wraparound=True, warn.unused=True, warn.unused_args=True,
# warn.unused_results=True

from bonndit.utilc.cython_helpers cimport sub_vectors, angle_deg, sum_c, set_zero_matrix, dist
import numpy as np
DTYPE = np.float64



cdef class Validator:
	def __cinit__(self, double[:,:,:] wm_mask, int[:] shape, double min_wm, double[:,:] inclusion, int r, double max_angle):
		self.min_wm = min_wm
		self.wm_mask = wm_mask
		self.shape = shape
		if r > 0:
			self.ROI = ROIValidator(inclusion, r)
		else:
			self.ROI = ROINotValidator(inclusion, r)
		if max_angle > 0:
			self.Curve = CurvatureValidator(max_angle)
		else:
			self.Curve = CurvatureNotValidator(max_angle)



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
	def __cinit__(self, max_angle):
		self.max_angle = max_angle
		self.angle = 0
		self.points = np.zeros([2,3])

	cdef bint curvature_checker(self, double[:,:] path, int path_len, double[:] features) nogil except *:
		return True

cdef class CurvatureValidator(CurvatureNotValidator):
	#def __cinit__(self, double max_angle):
#		super().__cinit__(max_angle)

	cdef bint curvature_checker(self, double[:,:] path, int path_len, double[:] features) nogil except *:
			"""
			Checks the angles between the current direction and the directions anlong the polygon. If a angle is to large returns True
			@param path: polygon to check
			@param path_len: length of the polygon
			@param features: save the angle between the current direction and the direction k points ago into the features.
			@return:
			"""
			cdef int l
			sub_vectors(self.points[1], path[path_len], path[path_len + 1])
			for l in range(path_len):
				sub_vectors(self.points[0], path[path_len - l], path[path_len - l + 1])
				if sum_c(self.points[1]) != 0 and sum_c(self.points[0]) != 0:
					self.angle = angle_deg(self.points[1], self.points[0])
				if self.angle > self.max_angle:
					return True
			else:
				features[0] = self.angle
				return False

cdef class ROINotValidator:
	def __cinit__(self, double[:,:] inclusion, int r):
		self.inclusion = inclusion
		self.inclusion_check = np.zeros(r)
		self.inclusion_num = r


	cdef void included(self, double[:] point) nogil except *:
		pass

	cdef bint included_checker(self) nogil except *:
		return False

cdef class ROIValidator(ROINotValidator):
	#def __cinit__(self, double[:,:] inclusion, int r ):
	#	super().__cinit__(inclusion, r)

	cdef void included(self, double[:] point) nogil except *:
		if sum_c(self.inclusion_check) == self.inclusion_num:
			return
		for i in range(self.inclusion.shape[0]):
			if dist(self.inclusion[i,1:], point) < 3:
				self.inclusion_check[int(self.inclusion[i,0])] = 1
				break

	cdef bint included_checker(self) nogil except *:
		return sum_c(self.inclusion_check) != self.inclusion_num
