#%%cython --annotate
#cython: language_level=3, boundscheck=False, wraparound=False, warn.unused=True, warn.unused_args=True



from .ItoW cimport Trafo
from .alignedDirection cimport Probabilities
from .kalman.model cimport AbstractModel
from .kalman.kalman cimport Kalman

cdef class Interpolation:
	cdef double[:,:,:,:,:] vector_field
	cdef double[:,:,:] cuboid
	cdef double[:,:]  floor_point
	cdef double[:,:] inv_trafo
	cdef double[:] point_index
	cdef double[:] point_world
	cdef double[:,:] best_dir
	cdef int[:,:,:,:] cache
	cdef double chosen_angle
	cdef Probabilities prob
	cdef int best_ind
	cdef double[:] next_dir, vector
	cdef void main_dir(self, double[:]) # nogil
	cdef void calc_cube(self, double[:]) # nogil
	cdef void nearest_neigh(self, double[:]) # nogil
	cdef void set_vector(self, int, int) # nogil
	cdef int interpolate(self, double[:], double[:], int) # nogil except *

cdef class FACT(Interpolation):
	cdef int interpolate(self, double[:], double[:], int) # nogil except *


cdef class Trilinear(Interpolation):
	cdef double[:,:] array, x_array, new_best_dir
	cdef int[:,:] not_check
	cdef double[:] point
	cdef double[:,:,:] dir
#	cdef int[:,:,:,:] cache
	cdef int[:] floor
	cdef int[:] permutation
	cdef void set_array(self, int, int, int) # nogil
	cpdef int interpolate(self, double[:], double[:], int) # nogil except *
	cdef void set_new_poss(self) # nogil except *
	cdef int kmeans(self, double[:]) # nogil except *
	cdef void permute(self, double[:]) # nogil except *

cdef class TrilinearFODF(Interpolation):
	cdef double[:,:,:,:] data
	cdef double[:] fodf
	cdef double[:] fodf1
	cdef double[:] empty
	cdef double sigma_1
	cdef double sigma_2
	cdef int inc
	cdef double[:] point_diff
	cdef double[:,:] trafo
	cdef double[:] dist
	cdef double[:] length
	cdef bint auto
	cdef double[:,:] best_dir_approx
	cdef double r
	cdef int rank
	cdef double[:,:] vlinear
	cdef int[:,:] neighbors
	cdef void trilinear(self, double[:] point) # nogil except *
	cdef void neigh(self, double[:] point) # nogil except *
	cdef int interpolate(self, double[:] point, double[:] old_dir, int r) # nogil except *

cdef class TrilinearFODFWatson(Interpolation):
	cdef double[:,:,:,:] data
	cdef double[:] fodf
	cdef double[:] fodf1
	cdef double[:] empty

	# cdef double[:] x_v
	# cdef double[:] dipy_v
	# cdef double[:] signal_v
	# cdef double[:] est_signal_v
	# cdef double[:,:,:] pysh_v
	# cdef double[:,:,:] rot_pysh_v
	# cdef double[:] angles_v
	# cdef double[:,:,:] dj_v
	# cdef double[:] loss_v

	cdef int lmax
	cdef int amount
	cdef double[:,:,:] dj
	cdef double[:,:] x_v
	cdef double[:,:] x_v2
	cdef double[:,:] signals
	cdef double[:] loss
	cdef double[:,:] angles_v
	cdef double[:,:] dipy_v
	cdef double[:,:,:,:] pysh_v
	cdef double[:,:,:,:] rot_pysh_v
	cdef double[:,:] est_signal
	cdef double[:] kappas
	cdef double[:] weights
	cdef double[:] kappa_range

	cdef double sigma_1
	cdef double sigma_2
	cdef int inc
	cdef double[:] point_diff
	cdef double[:,:] trafo
	cdef double[:] dist
	cdef double[:] length
	cdef bint auto
	cdef double[:,:] best_dir_approx
	cdef double r
	cdef int rank
	cdef double[:,:] vlinear
	cdef int[:,:] neighbors
	cdef void trilinear(self, double[:] point) # nogil except *
	cdef void neigh(self, double[:] point) # nogil except *
	cdef int interpolate(self, double[:] point, double[:] old_dir, int r) # nogil except *

cdef class UKF(Interpolation):
	cdef double[:] mean
	cdef double[:,:] P
	cdef double[:,:,:,:] data
	cdef double[:,:] mlinear
	cdef double[:] y
	cdef Kalman _kalman
	cdef AbstractModel _model
	cdef int interpolate(self, double[:], double[:], int) # nogil except *

cdef class UKFFodf(UKF):
	cdef int interpolate(self, double[:], double[:], int) # nogil except *

cdef class UKFMultiTensor(UKF):
	cdef int interpolate(self, double[:], double[:], int) # nogil except *









