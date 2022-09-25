#%%cython --annotate
#cython: language_level=3, boundscheck=False,
import cython
from libc.math cimport acos, pi, exp, fabs, cos, pow, tanh, round, sqrt
from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
from bonndit.utilc.cython_helpers cimport scalar, clip, mult_with_scalar, sum_c, norm, angle_deg
import numpy as np
from scipy.special import dawsn


###
# Given a set of vectors should be in shape [n,3] and a vector returns most aligned vector shape [3]

cdef class Probabilities:
	srand(time(NULL))



	def __cinit__(self, double expectation=0, double sigma=9):
		self.sigma = sigma
		self.expectation = expectation
		self.chosen_prob = 0
		self.probability = np.zeros((3,))
		self.angles = np.zeros((3,))
		self.best_fit = np.zeros((3))
		self.test_vectors = np.zeros((3, 3))
		self.chosen_angle = 0
		self.old_fa = 1


	cdef void aligned_direction(self, double[:,:] vectors, double[:] direction) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		#calculate angle between direction and possibilities. If angle is bigger than 90 use the opposite direction.
		cdef int i, n = vectors.shape[0]
		cdef double test_angle, min_angle = 180

		for i in range(n):
			if sum_c(direction) == 0:
				self.angles[i] = 0
				continue
			if norm(vectors[i]) != 0 and norm(vectors[i]) == norm(vectors[i]):
				test_angle = clip(scalar(direction, vectors[i])/(norm(direction)*(norm(vectors[i]))), -1,1)
				if test_angle >0 :
					self.angles[i] = acos(test_angle)/pi*180
					self.test_vectors[i] = vectors[i]
				elif test_angle <= 0:
					self.angles[i] = 180 - acos(test_angle)/pi*180
					mult_with_scalar(self.test_vectors[i], -1, vectors[i])
			else:
				self.angles[i] = 180
				mult_with_scalar(self.test_vectors[i], 0, vectors[i])

	cdef void random_choice(self, double[:] direction) : # nogil  except *:
		"""

		@param direction:
		@return:
		"""
		cdef double best_choice = rand() / RAND_MAX
	#	with gil:
	#		print(*self.probability)
		if sum_c(self.probability) != 0:
			mult_with_scalar(self.probability, 1/sum_c(self.probability), self.probability)


			if best_choice < self.probability[0]:
				mult_with_scalar(self.best_fit, 1, self.test_vectors[0])
				self.chosen_prob = self.probability[0]
				self.chosen_angle = self.angles[0]
			elif best_choice < self.probability[0] + self.probability[1]:
				mult_with_scalar(self.best_fit, 1, self.test_vectors[1])
				self.chosen_prob = self.probability[1]
				self.chosen_angle = self.angles[1]
			else:
				mult_with_scalar(self.best_fit, 1, self.test_vectors[2])
				self.chosen_prob = self.probability[2]
				self.chosen_angle = self.angles[2]
			self.old_fa = norm(self.best_fit)
		else:

			mult_with_scalar(self.best_fit, 0, self.test_vectors[2])
			self.chosen_angle = 0
			self.chosen_prob = 0
		#with gil:
		#	print(*self.probability, ' and I chose ', self.chosen_prob, ' where the angle are ', *self.angles, ' I chose ', self.chosen_angle)



	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil except *:
		pass

	cdef void calculate_watson_probabilities(self, double[:,:] vectors, double[:] kappas, double[:] weights, double[:] direction, double[:] point) : # nogil except *:
		self.calculate_probabilities(vectors, direction, point)

cdef class Gaussian(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i
		self.aligned_direction(vectors, direction)
		for i in range(3):
			self.probability[i] = exp(-1/2*((self.angles[i] - self.expectation)/self.sigma)**2)
		self.random_choice(direction)


cdef class Laplacian(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i
		self.aligned_direction(vectors, direction)
		for i in range(3):
			self.probability[i] = 1/2 * exp(- (fabs(self.angles[i] - self.expectation) /
			                                             self.sigma))
		self.random_choice(direction)


cdef class ScalarOld(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i
		cdef double s
		self.aligned_direction(vectors, direction)
		#with gil:
		#	print(*self.angles)
		for i in range(3):
			if sum_c(self.test_vectors[i]) == sum_c(self.test_vectors[i])  and pow(self.expectation/pow(2*pi,0.5)*self.angles[i]/180*pi,2) <= 1/2*pi:
			#	with gil:
			#		print('First angle ' , self.angles[i], pow(cos(pow(self.expectation/pow(2*pi,0.5)*self.angles[i]/180*pi,2)),self.sigma)*norm(self.test_vectors[i]))
				self.probability[i]=pow(cos(pow(self.expectation/pow(2*pi,0.5)*self.angles[i]/180*pi,2)),self.sigma)*exp(-pow(norm(self.test_vectors[i]) - self.old_fa,2)/0.01)


			else:
				self.probability[i] = 0
		self.random_choice(direction)


cdef class ScalarNew(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i
		cdef double s
		self.aligned_direction(vectors, direction)
		#with gil:
		#	print(*self.angles)
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]):
				self.probability[i] = pow(cos(self.angles[i]/180*pi),self.sigma)*norm(self.test_vectors[i])
			else:
				self.probability[i] = 0

		self.random_choice(direction)


cdef class Deterministic2(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i, min_index=0
		cdef double s, min_angle=0
		self.aligned_direction(vectors, direction)
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]) and sum_c(vectors[i])!=0:
				if self.angles[i] < min_angle or i==0:
					min_angle=self.angles[i]
					min_index=i
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]) and sum_c(vectors[i])!=0:
				if self.angles[i] < self.expectation or (self.angles[i] == min_angle and min_angle < self.sigma):
					self.probability[i] = 1
				else:
					self.probability[i] = 0
	#	self.probability[min_index] = 1
		self.random_choice(direction)
#		mult_with_scalar(self.best_fit, 1, self.test_vectors[min_index])
#		self.chosen_prob = 0
#		self.chosen_angle = self.angles[min_index]

cdef class Deterministic(Probabilities):
	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i, min_index=0
		cdef double s, min_angle=0
		self.aligned_direction(vectors, direction)
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]) and sum_c(vectors[i])!=0:
				if self.angles[i] < min_angle or i==0:
					min_angle=self.angles[i]
					min_index=i

		mult_with_scalar(self.best_fit, 1, self.test_vectors[min_index])
		self.chosen_prob = 0
		self.chosen_angle = self.angles[min_index]

cdef class Watson(Probabilities):
	cdef void watson_config(self, double[:,:,:,:] kappa_field) :#nogil  except *:
		self.kappa_field = kappa_field

	cdef double poly_kummer(self, double kappa) :#nogil  except *:
		return exp(kappa)/sqrt(kappa) * dawsn(sqrt(kappa))
		#return 3.52113645e+00 + kappa*(-3.13883527e+01 + kappa*(1.25022169e+02 + kappa*(-2.32549061e+02 + kappa*(2.47116824e+02 + kappa*(-1.66808904e+02 + kappa*(7.65540870e+01 + kappa*(-2.49902824e+01 + kappa*(5.98438651e+00 + kappa*(-1.07327464e+00 + kappa*(1.46034113e-01 + kappa*(-1.51685928e-02 + kappa*(1.20240009e-03 + kappa*(-7.21893666e-05 + kappa*(3.22742537e-06 + kappa*(-1.04170221e-07 + kappa*(2.29586418e-09 + kappa*(-3.09649360e-11 + kappa*(1.93236159e-13))))))))))))))))))

	cdef double poly_watson(self, double[:] x, double[:] mu, double kappa) :#nogil  except *:
		cdef double M = 4*pi*self.poly_kummer(kappa)
		return 1/M * exp(kappa * scalar(mu,x)**2)
	
	# rejection sampling from watson - watson_confidence_interval.ipynb
	cdef void mc_random_direction(self, double[:] direction, double[:] mu, double kappa) :#nogil  except *:
		cdef double max_val = self.poly_watson(mu, mu, kappa)
		cdef bint accept = False
		cdef double val, cutoff

		while not accept:
			direction[0] = (rand() / RAND_MAX) * 2 - 1
			direction[1] = (rand() / RAND_MAX) * 2 - 1
			direction[2] = (rand() / RAND_MAX) * 2 - 1
			mult_with_scalar(direction,1/norm(direction),direction)
			val = self.poly_watson(direction, mu, kappa)
			cutoff = (rand() / RAND_MAX) * max_val
			if val > cutoff:
				accept = True

	cdef void calculate_probabilities(self, double[:,:] vectors, double[:] direction, double[:] point) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i, min_index=0
		cdef double s, min_angle=0, norm_of_test

		cdef double[:] watson_point

		self.aligned_direction(vectors, direction)
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]) and sum_c(vectors[i])!=0:
				if self.angles[i] < min_angle or i==0:
					min_angle=self.angles[i]
					min_index=i

		# normalize length of selected peak direction
		norm_of_test = norm(self.test_vectors[min_index])
		if norm_of_test != 0:
			mult_with_scalar(self.test_vectors[min_index],1/norm_of_test,self.test_vectors[min_index])

		self.mc_random_direction(self.best_fit, 
								 self.test_vectors[min_index], 
								 #self.kappa_field[min_index, int(round(point[0])), int(round(point[1])), int(round(point[2]))])
								 max(10, min(50,self.kappa_field[min_index, int(round(point[0])), int(round(point[1])), int(round(point[2]))])))
								 #max(8, min(20,self.kappa_field[min_index, int(round(point[0])), int(round(point[1])), int(round(point[2]))]+8)))
		
		# flip direction if > 90:
		if scalar(self.best_fit, self.test_vectors[min_index]) < 0:
			mult_with_scalar(self.best_fit,-1.0,self.best_fit)

		#with gil:
			#print('Point', point[0], point[1], point[2], self.test)
			#print('Point', point[0], point[1], point[2], int(round(point[0])), int(round(point[1])), int(round(point[2])), self.kappa_field[0, int(round(point[0])), int(round(point[1])), int(round(point[2]))])
			#print(angle_deg(self.best_fit, self.test_vectors[min_index]), norm(self.test_vectors[min_index]), '   ', self.best_fit[0], self.best_fit[1], self.best_fit[2], '  ', self.test_vectors[min_index][0], self.test_vectors[min_index][1], self.test_vectors[min_index][2])

		#mult_with_scalar(self.best_fit, 1, self.test_vectors[min_index])

		# reset to original length
		mult_with_scalar(self.best_fit, norm_of_test, self.best_fit)

		#self.chosen_prob = ((max(8, min(20,self.kappa_field[min_index, int(round(point[0])), int(round(point[1])), int(round(point[2]))]+8))) - 8) / 12.
		self.chosen_prob = max(10, min(50,self.kappa_field[min_index, int(round(point[0])), int(round(point[1])), int(round(point[2]))]))
		self.chosen_angle = self.angles[min_index]

	cdef void calculate_watson_probabilities(self, double[:,:] vectors, double[:] kappas, double[:] weights, double[:] direction, double[:] point) : # nogil  except *:
		"""

		@param vectors:
		@param direction:
		@return:
		"""
		cdef int i, min_index=0
		cdef double s, min_angle=0, norm_of_test

		cdef double[:] watson_point

		self.aligned_direction(vectors, direction)
		for i in range(3):
			if sum_c(vectors[i]) == sum_c(vectors[i]) and sum_c(vectors[i])!=0:
				if self.angles[i] < min_angle or i==0:
					min_angle=self.angles[i]
					min_index=i

		# normalize length of selected peak direction
		norm_of_test = norm(self.test_vectors[min_index])
		mult_with_scalar(self.test_vectors[min_index],1/norm_of_test,self.test_vectors[min_index])

		self.mc_random_direction(self.best_fit, 
								 self.test_vectors[min_index], 
								 max(10, min(58,kappas[min_index])))
		
		# flip direction if > 90:
		if scalar(self.best_fit, self.test_vectors[min_index]) < 0:
			mult_with_scalar(self.best_fit,-1.0,self.best_fit)

		#with gil:
			#print('Point', point[0], point[1], point[2], self.test)
			#print('Point', point[0], point[1], point[2], int(round(point[0])), int(round(point[1])), int(round(point[2])), self.kappa_field[0, int(round(point[0])), int(round(point[1])), int(round(point[2]))])
			#print(angle_deg(self.best_fit, self.test_vectors[min_index]), norm(self.test_vectors[min_index]), '   ', self.best_fit[0], self.best_fit[1], self.best_fit[2], '  ', self.test_vectors[min_index][0], self.test_vectors[min_index][1], self.test_vectors[min_index][2])

		#mult_with_scalar(self.best_fit, 1, self.test_vectors[min_index])

		# reset to original length
		mult_with_scalar(self.best_fit, norm_of_test, self.best_fit)

		self.chosen_prob = max(10, min(58,kappas[min_index]))
		self.chosen_angle = self.angles[min_index]

