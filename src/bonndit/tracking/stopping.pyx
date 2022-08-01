#%%cython --annotate
#cython: language_level=3, boundscheck=False, wraparound=False, warn.unused=True, warn.unused_args=True,
# warn.unused_results=True
import os
from scipy.interpolate import RegularGridInterpolator

from bonndit.utilc.cython_helpers cimport sub_vectors, angle_deg, sum_c, set_zero_matrix, bigger, smaller, mult_with_scalar, norm
import numpy as np
from bonndit.tracking.ItoW cimport Trafo
DTYPE = np.float64


## TODO Mich um die Registrierung kümmern.
cdef class Validator:
	def __cinit__(self, double[:,:,:] wm_mask, int[:] shape, double min_wm, inclusion, exclusion,  double
		max_angle, Trafo trafo, double step_width, **kwargs):
		#self.min_wm = min_wm
		#self.wm_mask = wm_mask
		self.shape = shape
		if 'act' in kwargs:
			self.WM = ACT(kwargs['act'])
		else:
			self.WM = WMChecker(wm_mask, min_wm)
		if isinstance(inclusion, np.ndarray):
			self.ROIIn = ROIInValidator(inclusion)
		else:
			self.ROIIn = ROIInNotValidator(np.zeros((3,3)))
		if isinstance(exclusion, np.ndarray):
			self.ROIEx = ROIExValidator(exclusion)
		else:
			self.ROIEx = ROIExNotValidator(np.zeros((3,3)))
		if max_angle > 0:
			self.Curve = CurvatureValidator(max_angle, trafo, step_width)
		else:
			self.Curve = CurvatureNotValidator(max_angle, trafo, step_width)



	#cdef bint wm_checker(self, double[:] point) : # nogil except *:
	#	""" Checks if the wm density is at a given point below a threshold.
	#	@param point: 3 dimensional point
	#	"""
	#	if self.wm_mask[int(point[0]), int(point[1]), int(point[2])] < self.min_wm:
	#		return True
	#	else:
	#		return False

	cdef bint index_checker(self, double[:] point) : # nogil except *:
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



	cdef bint next_point_checker(self, double[:] point) : # nogil except *:
		"""
		Check if a given direction is valid e.g. not zero and not infinity
		@param point: given direction
		@return:
		"""
		if sum_c(point) == 0 or sum_c(point) != sum_c(point):
			return True
		else:
			return False

	cdef void set_path_zero(self, double[:,:] path, double[:,:] features) : # nogil except *:
		set_zero_matrix(path)
		set_zero_matrix(features)

cdef class WMChecker:
	cdef __cinit__(self, wm_mask, min_wm):
		self.min_wm = min_wm
		self.wm_mask = wm_mask

	cdef void reset(self):
		pass

	cdef bint sgm_checker(self, double[:] point):
		return 0

	cdef bint wm_checker(self, double[:] point) : # nogil except *:
			""" Checks if the wm density is at a given point below a threshold.
			@param point: 3 dimensional point
			"""
			if self.wm_mask[int(point[0]), int(point[1]), int(point[2])] < self.min_wm:
				return 0
			else:
				return -1

cdef class ACT:
	"""
	Format of five_tt:
    0: Cortical grey matter
    1: Sub-cortical grey matter
    2: White matter
    3: CSF
    4: Pathological tissue
	"""
	cdef __cinit__(self, five_tt):
		x = np.linspace(0, five_tt.shape[1], five_tt.shape[1])
		y = np.linspace(0, five_tt.shape[2], five_tt.shape[2])
		z = np.linspace(0, five_tt.shape[3], five_tt.shape[3])
		self.entered_sgm = 0
		self.cgm = RegularGridInterpolator((x,y,z), five_tt[0])
		self.sgm = RegularGridInterpolator((x, y, z), five_tt[1])
		self.wm = RegularGridInterpolator((x, y, z), five_tt[2])
		self.csf = RegularGridInterpolator((x, y, z), five_tt[3])

	cdef void reset(self):
		self.entered_sgm = 0

	cdef bint wm_checker(self, double[:] point):
		cgm = self.cgm(point)
		csf = self.csf(point)
		sgm = self.sgm(point)
		wm = self.sgm(point)
		#check case 6:
		if wm > 0.5:
			return -1
		if self.entered_sgm:
			if sgm<0.5:
				return 1
		# ACT cases 1..3
		if cgm > 0.5:
			return 1
		if csf > 0.5:
			return 2
		if cgm + csf + sgm + wm < 0.3:
			return 1
		#case 6
		if sgm>0.5:
			self.entered_sgm = 1




	cdef bint sgm_checker(self, double[:] point):
		sgm = self.sgm(point)
		if sgm>0.5:
			return 1
		else:
			return 0


cdef class CurvatureNotValidator:
	def __cinit__(self, max_angle, trafo, double step_width):
		self.max_angle = max_angle
		self.angle = 0
		self.step_width = step_width
		self.points = np.zeros([5,3])
		self.trafo = trafo

	cdef bint curvature_checker(self, double[:,:] path,  double[:] features) : # nogil except *:
		return False

cdef class CurvatureValidator(CurvatureNotValidator):
	#def __cinit__(self, double max_angle):
#		super().__cinit__(max_angle)

	cdef bint curvature_checker(self, double[:,:] path,  double[:] features) : # nogil except *:
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
			sub_vectors(self.points[1], path[k-1], path[k])
			length += norm(self.points[1])
			#with gil: print(k, length, l)
			while k >= 2 and length < 30 and l < k:
				l += 1
				sub_vectors(self.points[0], path[k-l], path[k-l+1])
				length += norm(self.points[0])
				self.angle = angle_deg(self.points[1], self.points[0])
			#	with gil:
			#		print("test", self.angle)
				if self.angle > self.max_angle:
					return True
			else:
				features[0] = self.angle
				return False

cdef class ROIInNotValidator:
	def __cinit__(self, double[:,:] inclusion):
		self.inclusion = np.zeros([3,3])
		self.inclusion_num = 0
		self.inclusion_check = np.zeros(1)


	cdef int included(self, double[:] point) : # nogil except *:
		return 0

	cdef bint included_checker(self) : # nogil except *:
		return False

cdef class ROIInValidator(ROIInNotValidator):
	def __cinit__(self, double[:,:] inclusion):
		self.inclusion = inclusion[:,:3]
		self.inclusion_num = inclusion.shape[0]//2
		self.inclusion_check = np.zeros(inclusion.shape[0]//2)

	cdef int included(self, double[:] point) : # nogil except *:
		cdef int i
		if sum_c(self.inclusion_check) == self.inclusion_num:
			return 0

		for i in range(self.inclusion_num):
			if not bigger(point, self.inclusion[2*i]):
				continue
			if not smaller(point, self.inclusion[2*i + 1]):
				continue
			self.inclusion_check[i] = 1
			return i + 1
		return -1

	cpdef int included_p(self, double[:] point) except *:
		return self.included(point)

	cpdef void reset_p(self) except *:
		self.inclusion_check = np.zeros(self.inclusion_num)


	cdef bint included_checker(self) : # nogil except *:
		return sum_c(self.inclusion_check) != self.inclusion_num


cdef class ROIExNotValidator:
	def __cinit__(self, double[:,:] exclusion):
		self.exclusion_cube = np.zeros([3,3])
		self.exclusion_num = 0

	cdef bint excluded(self, double[:] point) : # nogil except *:
		return False


cdef class ROIExValidator(ROIExNotValidator):
	def __cinit__(self, double[:,:] exclusion):
		self.exclusion_cube = exclusion
		self.exclusion_num = len(exclusion.shape[0])

	cdef bint excluded(self, double[:] point) : # nogil except *:
		cdef int i
		for i in range(self.exclusion_num):
			if not bigger(point, self.exclusion_cube[2*i]):
				continue
			if not smaller(point, self.exclusion_cube[2*i + 1]):
				continue
			return True
		return False




