from abstract_adapt_operation import *

class WeightedRatio(DimensionAdaptivity):

	def __init__(self, dim, tol, init_multiindex, max_level, w):
		
		self._dim 				= dim
		self._tol 				= tol
		self._init_multiindex 	= init_multiindex
		self._max_level 		= max_level
		self.__w 				= w 

		self._approx 	= 0.
		self._eta 		= 0.

		self._O 				= OrderedDict()
		self._A 				= OrderedDict() 
		self._local_error 		= OrderedDict()
		self._key_O 			= 0
		self._key_A 			= 0
		self._key_local_error 	= 0

		self._multiindex_set = []

		self._init_delta 		= 0.
		self._init_no_points 	= 0

		self._stop_adaption = False

	def _get_norm_delta(self, delta):

		return np.abs(delta)

	def _get_local_error_idicator(self, init_delta, curr_delta, init_no_points, curr_no_points):
		
		local_error = max(self.__w * self._get_norm_delta(curr_delta)/self._get_norm_delta(init_delta), \
											(1. - self.__w) * init_no_points/curr_no_points)

		return local_error

	def init_adaption(self, init_delta, init_no_points):

		self._key_O 								= -1
		self._key_A 								= 0
		self._key_local_error 						= 0
		self._A[self._key_A] 						= self._init_multiindex
		local_error_indicator 						= self._get_local_error_idicator(init_delta, init_delta, \
																init_no_points, init_no_points)
		self._local_error[self._key_local_error] 	= local_error_indicator

		self._init_delta 		= init_delta
		self._init_no_points 	= init_no_points

		self._approx = init_delta
		self._eta 	 = local_error_indicator

		self._multiindex_set.append(self._init_multiindex)

	def do_one_adaption_step_preproc(self):

		local_multiindices = []

		max_index 	= self._select_highest_priority_index()
		max_i 		= self._A[max_index]
		
		self._key_O				+= 1
		self._O[self._key_O] 	= max_i
		self._eta 				-= self._local_error[max_index]

		del self._A[max_index]
		del self._local_error[max_index]

		neighbors_i = Multiindex(self._dim).get_successors(max_i)
		for neighbor in neighbors_i:
			if self._is_O_admissible(neighbor):
				local_multiindices.append(neighbor)

				self._key_A 		+= 1
				self._A[self._key_A] = neighbor

				self._multiindex_set.append(neighbor)

		local_multiindices = np.array(local_multiindices, dtype=int)

		return local_multiindices

	def do_one_adaption_step_postproc(self, deltas, no_points):

		for no_points_level, delta_level in zip(no_points, deltas):
			self._key_local_error 					+= 1
			local_error_indicator					 = self._get_local_error_idicator(self._init_delta, delta_level, \
																	self._init_no_points, no_points_level)
			self._local_error[self._key_local_error] = local_error_indicator

			self._approx 	+= delta_level
			self._eta 		+= local_error_indicator

	def check_termination_criterion(self):

		max_level = np.max(self._multiindex_set)
		if len(self._A.values()) == 0 or self._eta <= self._tol or max_level >= self._max_level:
			self._stop_adaption = True 