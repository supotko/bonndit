#cython: language_level=3, boundscheck=False, wraparound=False, warn.unused=True, warn.unused_args=True
# warn.unused_results=True
import sys
import nrrd
sys.path.append('.')
from .alignedDirection cimport  Gaussian, Laplacian, ScalarOld, ScalarNew, Probabilities, Deterministic,Deterministic2, Watson
from .ItoW cimport Trafo
from .stopping cimport Validator
from .integration cimport  Euler, Integration, EulerUKF, RungeKutta
from .interpolation cimport  FACT, Trilinear, Interpolation, UKFFodf, UKFMultiTensor, TrilinearFODF, TrilinearFODFWatson
from bonndit.utilc.cython_helpers cimport mult_with_scalar, sum_c, sum_c_int, set_zero_vector, sub_vectors, \
	angle_deg, norm
import numpy as np
from tqdm import tqdm
from bonndit.utilc.myext cimport mw_openmp_single_o4c, mw_openmp_multc, mw_openmp_mult_o4c

ctypedef struct possible_features:
	int chosen_angle
	int seedpoint
	int prob_chosen
	int prob_others_0
	int prob_others_1
	int prob_others_2
	int fa
	int len


cdef void tracking(double[:,:,:,:] paths, double[:] seed,
                   int seed_shape, Interpolation interpolate,
              Integration integrate, Trafo trafo, Validator validator, int max_track_length, int save_steps,
	                   int samples, double[:,:,:,:] features, possible_features features_save, int minlen) : # nogil except *:
	"""
        Initializes the tracking for one seed in both directions.
	@param paths:
	@param old_dir:
	@param seed:
	@param seed_shape:
	@param interpolate:
	@param integrate:
	@param trafo:
	@param validator:
	@param max_track_length:
	@param samples:
	@param features:
	"""
	cdef int k=0, j, l, m, u

	for j in range(samples):
		k=0
		while True:
			k+=1
			# set zero inclusion check
			set_zero_vector(validator.ROIIn.inclusion_check)
			if seed_shape == 3:
				interpolate.main_dir(paths[j, 0, 0])
				integrate.old_dir = interpolate.next_dir
			else:
				integrate.old_dir = seed[3:]
			status1, m = forward_tracking(paths[j,:,0, :], interpolate, integrate, trafo, validator, max_track_length, save_steps,
			                 features[j,:,0, :], features_save,)

			if seed_shape == 3:
				interpolate.main_dir(paths[j, 0, 1])
				mult_with_scalar(integrate.old_dir, -1.0 ,interpolate.next_dir)
			else:
				mult_with_scalar(integrate.old_dir, -1.0 ,seed[3:])
			status2, l = forward_tracking(paths[j,:,1,:], interpolate, integrate, trafo, validator, max_track_length, save_steps,
			                 features[j,:,1, :], features_save)
			# if not found bot regions of interest delete path.
			if validator.ROIIn.included_checker() or not status1 or not status2 or (l+m)*save_steps*integrate.stepsize < minlen:
				validator.set_path_zero(paths[j,:,1,:], features[j,:,1, :])
				validator.set_path_zero(paths[j, :, 0, :], features[j, :, 0, :])
				for u in range(3):
					paths[j, 0, 0,u] = seed[u]
					paths[j, 0, 1,u] = seed[u]
				if features_save.seedpoint >= 0:
					features[j, 0, 0, features_save.seedpoint] = 1
					features[j, 0, 1, features_save.seedpoint] = 1
			else:
				break
			if k==2:
				break

cdef forward_tracking(double[:,:] paths,  Interpolation interpolate,
                       Integration integrate, Trafo trafo, Validator validator, int max_track_length, int save_steps, double[:,:] features, possible_features feature_save) : # nogil except *:
	"""
        This function do the tracking into one direction.
	@param paths: empty path array to save the streamline points.
    @param interpolate: Interpolation object
    @param integrate: Integration object
    @param trafo: Trafo object
    @param validator: Validator object
	@param max_track_length:
    @param features: empty feature array. To save informations to the streamline
	"""
	# check wm volume
	cdef int k, con
	# thousand is max length for pathway
	interpolate.prob.old_fa = 1
	validator.WM.reset()
	for k in range((max_track_length-1)*save_steps):
		# validate index and wm density.

		if validator.index_checker(paths[(k-1)//save_steps + 1]):
			set_zero_vector(paths[(k-1)//save_steps + 1])
			break
		# check if neigh is wm.
		con = validator.WM.wm_checker(paths[(k - 1) // save_steps + 1])
		if con == 0:
			# If not wm. Set all zero until wm
			for l in range((k - 1) // save_steps + 1):
				con = validator.WM.wm_checker_ex(paths[(k - 1) // save_steps + 1 - save_steps * l])
				if con == 0:
					set_zero_vector(paths[(k - 1) // save_steps + 1 - save_steps * l])
				else:
					break
			break
		elif con > 0:
			return False, k

		# find matching directions
		if sum_c(integrate.old_dir) == 0:
			set_zero_vector(paths[(k-1)//save_steps + 1])
			set_zero_vector(features[(k - 1) // save_steps + 1])

			break

		if interpolate.interpolate(paths[(k-1)//save_steps + 1], integrate.old_dir, (k-1)//save_steps + 1) != 0:
			break

		# Check next step is valid. If it is: Integrate. else break
		if validator.next_point_checker(interpolate.next_dir):
			set_zero_vector(paths[(k-1)//save_steps + 1])
			break

		if integrate.integrate(interpolate.next_dir, paths[(k-1)//save_steps + 1])!= 0:
			break

		if sum_c(integrate.next_point) == 0:
			break

		# update old dir
		paths[k//save_steps + 1] = integrate.next_point
		# check if next dir is near region of interest:

		validator.ROIIn.included(paths[k//save_steps])
		if validator.ROIEx.excluded(paths[k//save_steps]):
			return False, k
		if feature_save.chosen_angle >= 0:
			features[k//save_steps,feature_save.chosen_angle] = interpolate.prob.chosen_angle
		if feature_save.prob_chosen >= 0:
			features[k//save_steps,feature_save.prob_chosen] = interpolate.prob.chosen_prob
		if feature_save.prob_others_0 >= 0:
			features[k//save_steps,feature_save.prob_others_0] = interpolate.prob.probability[0]
			features[k//save_steps,feature_save.prob_others_1] = interpolate.prob.probability[1]
			features[k//save_steps,feature_save.prob_others_2] = interpolate.prob.probability[2]
		if feature_save.fa >= 0:
			features[k//save_steps,feature_save.fa] = interpolate.prob.old_fa
		# Check curvature between current point and point 30mm ago
		if validator.Curve.curvature_checker(paths[:k//save_steps], features[k//save_steps:k//save_steps + 1,1]):

			if not validator.WM.sgm_checker(trafo.point_wtoi):
				return False, k
		#integrate.old_dir = interpolate.next_dir

	#if k == 0:
	#	trafo.itow(paths[k//save_steps])
	#	paths[k//save_steps] = trafo.point_itow
	return True, k

cpdef tracking_all(vector_field, wm_mask, seeds, tracking_parameters, postprocessing, ukf_parameters, trilinear_parameters, logging, saving):
	"""
	@param vector_field: Array (4,3,x,y,z)
		Where the first dimension contains the length and direction, the second
		contains the directions.
	@param meta: Dictionary. Keys: space directions and space origin
		meta data of multivectorfield
	@param wm_mask: Array (x,y,z)
		WM Mask to check wm density
	@param seeds: Array (3,r) or (6,r)
		seedpoint array. If first first axis has length 3 only seed point is given.
		If 6 also the initial direcition is given.
	@param integration: String
		Only Euler integration is available.
	@param interpolation: String
		FACT und Trilinear are available
	@param prob: String
		Gaussian, Laplacian, ScalarOld and ScalarNew are implemented
	@param stepsize: double
		Stepsize for Euler integration.
	@param variance: double
		Variance of Probabilistic Selection Method.
	@param samples: int
		Count of samples per seed.
	@param max_track_length: int
		Maximal length of each track.
	@param wmmin: double
		Minimal White Matter density
	@param expectation: double
		Expectation of Probabilistic selection method.
	@return:
	"""
	cdef Interpolation interpolate
	cdef Integration integrate
	cdef Trafo trafo
	cdef Probabilities directionGetter
	cdef Validator validator

	cdef Watson temp_watson
	cdef double[:,:,:,:] kappa_field


	# cdef double[:] dipy_v = np.zeros((15,))
	# cdef double[:] est_signal_v = np.zeros((15,))
	# cdef double[:] x_v = np.array([ 		2.45448036,  3.81193276,  0.21720971, -0.01643096,  4.84025163,
    # 								 		3.43986788,  1.38213121,  1.40086381,  3.80655244,  4.42389366,
    # 								 		1.24552913, -2.84815265]) #nonzero_csd_init[:amount].copy()
	# cdef double[:] signal_v = np.array([	0.0667792 , -0.00764293, -0.00308147, -0.00726303, -0.00556099,
    # 								 		-0.00397964, -0.01252561, -0.01570279,  0.00730777, -0.03085464,
    # 								 		0.05140559, -0.00341419,  0.01492144,  0.01563105,  0.0043892 ,
    # 								 		-0.0071501 , -0.0063297 ,  0.00019099,  0.00813005,  0.00532886,
    # 								 		-0.01496287,  0.00619936, -0.01115468, -0.00603983, -0.01041317,
    # 								 		0.00047331, -0.00409659, -0.0191022 ])
	# cdef double[:,:,:] pysh_v = np.zeros((2,5,5))
	# cdef double[:,:,:] rot_pysh_v = np.zeros((2,5,5))
	# cdef double[:] angles_v = np.zeros((3,))
	# cdef double[:,:,:] dj_v = np.zeros((5,5,5))
	# cdef double[:] loss_v = np.zeros((1,))
	# print("x before",x_v[0],x_v[1],x_v[2],x_v[3],"fodf",signal_v[1],signal_v[2],signal_v[3])
	# mw_openmp_single_o4c(x_v, signal_v, est_signal_v, dipy_v, pysh_v, rot_pysh_v, angles_v, loss_v, 3)
	# print("x after",x_v[0],x_v[1],x_v[2],x_v[3],"loss",loss_v[0])

	num_of_dir = 3

	# cdef int amount = 1
	# cdef double[:,:,:] dj = np.zeros((5,5,5))
	# cdef double[:,:] x_v = np.array([[ 2.45448036,  3.81193276,  0.21720971, -0.01643096,  4.84025163,3.43986788,  1.38213121,  1.40086381,  3.80655244,  4.42389366,1.24552913, -2.84815265]]) #nonzero_csd_init[:amount].copy()
	# cdef double[:,:] signals = np.array([[ 0.0667792 , -0.00764293, -0.00308147, -0.00726303, -0.00556099,-0.00397964, -0.01252561, -0.01570279,  0.00730777, -0.03085464, 0.05140559, -0.00341419,  0.01492144,  0.01563105,  0.0043892 ,-0.0071501 , -0.0063297 ,  0.00019099,  0.00813005,  0.00532886,-0.01496287,  0.00619936, -0.01115468, -0.00603983, -0.01041317, 0.00047331, -0.00409659, -0.0191022 ]])
	# cdef double[:] loss = np.zeros(amount)
	# cdef double[:,:] angles_v = np.zeros((amount, 3))
	# cdef double[:,:] dipy_v = np.zeros((amount, 28))
	# cdef double[:,:,:,:] pysh_v = np.zeros((amount, 2, 7, 7))
	# cdef double[:,:,:,:] rot_pysh_v = np.zeros_like(pysh_v)
	# cdef double[:,:] est_signal = np.zeros((amount, 28))

	cdef int amount = 1
	cdef double[:,:,:] dj = np.zeros((5,5,5))
	cdef double[:,:] x_v = np.array([[ 2.45448036,  3.81193276,  0.21720971, -0.01643096,  4.84025163,3.43986788,  1.38213121,  1.40086381,  3.80655244,  4.42389366,1.24552913, -2.84815265]]) #nonzero_csd_init[:amount].copy()
	cdef double[:,:] signals = np.array([[ 0.0667792 , -0.00764293, -0.00308147, -0.00726303, -0.00556099,-0.00397964, -0.01252561, -0.01570279,  0.00730777, -0.03085464, 0.05140559, -0.00341419,  0.01492144,  0.01563105,  0.0043892 ]])
	cdef double[:] loss = np.zeros(amount)
	cdef double[:,:] angles_v = np.zeros((amount, 3))
	cdef double[:,:] dipy_v = np.zeros((amount, 15))
	cdef double[:,:,:,:] pysh_v = np.zeros((amount, 2, 5, 5))
	cdef double[:,:,:,:] rot_pysh_v = np.zeros_like(pysh_v)
	cdef double[:,:] est_signal = np.zeros((amount, 15))

	print("x before",x_v[0,0],x_v[0,1],x_v[0,2],x_v[0,3])
	mw_openmp_mult_o4c(x_v, signals, est_signal, dipy_v, pysh_v, rot_pysh_v, angles_v, dj, loss, amount, 4, 3)
	print("x after",x_v[0,0],x_v[0,1],x_v[0,2],x_v[0,3],"loss",loss[0])

	# modify vector_field if Watson to split into directions and kappa
	if tracking_parameters['prob'] == "Watson":
		kappa_field = vector_field[0,:,:,:,:]
		vector_field = vector_field[1:,:,:,:,:]

	#select appropriate model
	if tracking_parameters['prob'] == "Gaussian":
		directionGetter = Gaussian(0, tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "Laplacian":
		directionGetter = Laplacian(0, tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "ScalarOld":
		directionGetter = ScalarOld(tracking_parameters['expectation'], tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "ScalarNew":
		directionGetter = ScalarNew(tracking_parameters['expectation'], tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "Deterministic":
		directionGetter = Deterministic(tracking_parameters['expectation'], tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "Deterministic2":
		directionGetter = Deterministic2(tracking_parameters['expectation'], tracking_parameters['variance'])
	elif tracking_parameters['prob'] == "Watson":
		temp_watson = Watson(tracking_parameters['expectation'], tracking_parameters['variance'])
		temp_watson.watson_config(kappa_field)
		directionGetter = temp_watson
	else:
		logging.error('Gaussian or Laplacian or Scalar are available so far. ')
		return 0

	trafo = <Trafo> Trafo(np.float64(tracking_parameters['space directions']), np.float64(tracking_parameters['space origin']))
	trafo_matrix = np.zeros((4,4))
	trafo_matrix[:3,:3] = tracking_parameters['space directions']
	trafo_matrix[:3,3] = tracking_parameters['space origin']
	trafo_matrix[3,3] = 1
	validator = Validator(np.array(wm_mask.shape, dtype=np.intc), postprocessing['inclusion'], postprocessing['exclusion'], trafo,  **tracking_parameters)



	cdef int[:] dim = np.array(vector_field.shape, dtype=np.int32)

	if tracking_parameters['ukf'] == "MultiTensor":
		interpolate = UKFMultiTensor(vector_field, dim[2:5], directionGetter, **ukf_parameters)
	elif tracking_parameters['ukf'] == "LowRank":
		interpolate = UKFFodf(vector_field, dim[2:5], directionGetter, **ukf_parameters)
	elif tracking_parameters['interpolation'] == "FACT":
		interpolate = FACT(vector_field, dim[2:5], directionGetter, **tracking_parameters)
	elif tracking_parameters['interpolation'] == "Trilinear":
		interpolate = Trilinear(vector_field, dim[2:5], directionGetter, **tracking_parameters)
	elif tracking_parameters['interpolation'] == "TrilinearFODF":
		interpolate = TrilinearFODF(vector_field, dim[2:5], directionGetter, **trilinear_parameters)
	elif tracking_parameters['interpolation'] == "TrilinearFODFWatson":
		interpolate = TrilinearFODFWatson(vector_field, dim[2:5], directionGetter, **trilinear_parameters)
	else:
		logging.error('FACT, Triliniear or UKF for MultiTensor and low rank approximation are available so far.')
		return 0

	if tracking_parameters['ukf'] == "MultiTensor":
		integrate = EulerUKF(tracking_parameters['space directions'], tracking_parameters['space origin'], trafo,
							 float(tracking_parameters['stepsize']))
	elif tracking_parameters['integration'] == "Euler":
		integrate = Euler(tracking_parameters['space directions'], tracking_parameters['space origin'], trafo, float(tracking_parameters['stepsize']))
	elif tracking_parameters['integration'] == "RungeKutta":
		integrate = RungeKutta(tracking_parameters['space directions'], tracking_parameters['space origin'], trafo, float(tracking_parameters['stepsize']), **{'interpolate': interpolate})
	else:
		logging.error('Only Euler is available so far. Hence set Euler as argument.')
		return 0

	cdef int i, j, k, l, m = seeds.shape[0]
	# Array to save Polygons
	cdef double[:,:,:,:,:] paths = np.zeros((1 if saving['file'] else m, tracking_parameters['samples'], tracking_parameters['max_track_length'], 2, 3),dtype=np.float64)
	# Array to save features belonging to polygons
	cdef double[:,:,:,:,:] features = np.zeros((1 if saving['file'] else m, tracking_parameters['samples'], tracking_parameters['max_track_length'], 2, saving['features']['len']),dtype=np.float64)
	# loop through all seeds.
	tracks = []
	tracks_len = []
	k = 0

	for i in tqdm(range(m), disable=not tracking_parameters['verbose']):
		#k = 0 if saving['file'] else k+=1
		if saving['file']:
			k = 0
		else: k = k +1
		#Convert seedpoint
		#trafo.wtoi(seeds[i][:3])
		for j in range(tracking_parameters['samples']):
			validator.set_path_zero(paths[k,j, :, 1, :], features[k,j, :, 1, :])
			validator.set_path_zero(paths[k,j, :, 0, :], features[k,j, :, 0, :])

			for l in range(3):
				paths[k,j, 0, 0,l] = seeds[i][l]
				paths[k,j, 0, 1,l] = seeds[i][l]
		#	if "Deterministic" in tracking_parameters['prob'] or tracking_parameters['ukf'] == "LowRank":
			for l in range(3):
				paths[k,j, 0, 0,l] +=  np.random.normal(0,1)
				paths[k,j, 0, 1,l] = paths[k,j, 0, 0,l]

		if saving['features']['seedpoint'] >= 0:
			features[k,:, 0, 0, saving['features']['seedpoint']] = 1
			features[k,:, 0, 1, saving['features']['seedpoint']] = 1
	
		# if interpolation is set to TrilinearFODFWatson create new interpolation object
		#if tracking_parameters['interpolation'] == "TrilinearFODFWatson":
		#	interpolate = TrilinearFODFWatson(vector_field, dim[2:5], directionGetter, **trilinear_parameters)

	#	print("1", np.asarray(features[k,j,:,0]))
	#	print("1", np.asarray(features[k,j,:,1]))
		#Do the tracking for this seed with the direction
		tracking(paths[k], seeds[i], seeds[i].shape[0], interpolate, integrate, trafo, validator, tracking_parameters['max_track_length'], tracking_parameters['sw_save'], tracking_parameters['samples'], features[k], saving['features'], tracking_parameters['min_len'])
		# delete all zero arrays.

		for j in range(tracking_parameters['samples']):
			path = np.concatenate((np.asarray(paths[k,j]),np.asarray(features[k,j])), axis=-1)
		# seedpoint would be twice if first index is not skipped.
			path = np.concatenate((path[1:,0][::-1], path[:,1]))
			try:
				to_exclude = np.all(path[:,:3] == 0, axis=1)

				path = path[~to_exclude]
				if path.size == 0:
					continue
			#	print(path)
				if path.shape[0]>5:
					# Work on disk or ram. Ram might be faster but for large files disk is preferable.
					if saving['file']:
						with open(saving['file'] + 'len', 'a') as f:
							f.write(str(path.shape[0]) +'\n')
						with open(saving['file'], 'a') as f:
							for l in range(path.shape[0]):
								f.write(' '.join(map(str, path[l])) + "\n")
					else:
						tracks_len.append(path.shape[0])
						tracks += [tuple(x) for x in path]
			except:
				pass

	return tracks, tracks_len


