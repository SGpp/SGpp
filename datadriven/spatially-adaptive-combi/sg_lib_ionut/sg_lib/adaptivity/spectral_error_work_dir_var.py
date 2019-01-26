from abstract_adapt_operation import *

class SpectralErrorWorkDirVar(DimensionAdaptivity):

	def __init__(self, dim, tol, tol_dims, init_multiindex, max_level, level_to_nodes):
		
		self._dim 				= dim
		self._tol 				= tol
		self._init_multiindex 	= init_multiindex
		self._max_level 		= max_level
		self._tol_dims 			= tol_dims
		self._level_to_nodes 	= level_to_nodes

		self._eta = 0.

		self._O 				= OrderedDict()
		self._A 				= OrderedDict() 
		self._local_error 		= OrderedDict()
		self._key_O 			= 0
		self._key_A 			= 0
		self._key_local_error 	= 0

		self._multiindex_set 	= []
		self._init_no_points 	= 0

		self._levels_converged_dims 	= None
		self._levels_converged_dims_all = None
		self._partial_var_active_set 	= None

		self._stop_adaption = False

		self._local_basis_global 	= None
		self._local_basis_local 	= OrderedDict()

	def __is_less_or_equal_then(self, multiindex_a, multiindex_b):

		assert len(multiindex_a) == len(multiindex_b)

		is_less_or_equal_then = True

		for i in xrange(len(multiindex_a)):
			if multiindex_a[i] > multiindex_b[i]:
				is_less_or_equal_then = False
				break

		return is_less_or_equal_then

	def _get_norm_delta(self, delta_coeff):

		norm = 0.

		if not isinstance(delta_coeff, np.ndarray) == 1:
			norm = np.asscalar(delta_coeff)**2
		else:
			norm = np.sum([c**2 for c in delta_coeff])

		norm = np.sqrt(norm)

		return norm

	def _get_local_error_idicator(self, delta_coeff, no_points):
		
		local_error = self._get_norm_delta(delta_coeff)/no_points

		return local_error

	@property
	def partial_var_active_set(self):
		
		return self._partial_var_active_set

	@property
	def levels_converged_dims(self):
		
		return self._levels_converged_dims

	@property
	def levels_converged_dims_all(self):
		
		return self._levels_converged_dims_all

	def init_adaption(self, init_coeff, init_no_points):

		self._key_O 								= -1
		self._key_A 								= 0
		self._key_local_error 						= 0
		self._A[self._key_A] 						= self._init_multiindex
		local_error_indicator 						= self._get_local_error_idicator(init_coeff, init_no_points)
		self._local_error[self._key_local_error] 	= local_error_indicator
		self._eta 	 								= local_error_indicator

		self._init_no_points = init_no_points
		self._multiindex_set.append(self._init_multiindex)

		self._levels_converged_dims 	= self._max_level*np.ones(self._dim, dtype=int)
		self._levels_converged_dims_all = []
		self._partial_var_active_set 	= np.ones(self._dim)

		self._local_basis_local[repr(self._init_multiindex)] 	= self._get_local_hierarchical_basis(self._init_multiindex)
		self._local_basis_global 								= self._get_local_hierarchical_basis(self._init_multiindex)

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
			if self._is_O_admissible(neighbor) and self.__is_less_or_equal_then(neighbor, self._levels_converged_dims):

				local_multiindices.append(neighbor)

				self._key_A 		+= 1
				self._A[self._key_A] = neighbor

				self._multiindex_set.append(neighbor)

				local_basis_neighbor = np.array([self._get_no_1D_grid_points(n) - 1 for n in neighbor], dtype=int)
				self._update_local_basis(neighbor.tolist(), local_basis_neighbor)

		local_multiindices = np.array(local_multiindices, dtype=int)

		return local_multiindices

	def do_one_adaption_step_postproc(self, curr_coeffs, no_points, spectral_op_obj):

		for no_points_level, delta_coeff in zip(no_points, curr_coeffs):
			self._key_local_error 					+= 1
			local_error_indicator					 = self._get_local_error_idicator(delta_coeff, no_points_level)
			self._local_error[self._key_local_error] = local_error_indicator

			self._eta += local_error_indicator

		self._partial_var_active_set 	= spectral_op_obj.get_directional_var_active_set(self._multiindex_set, self._A.values())
		curr_active_set 				= np.array(self._A.values())

		if len(curr_active_set) >= 1:
			curr_max_levels = np.amax(self._multiindex_set, axis=0)

			if len(self._levels_converged_dims_all) == 0:
				for d in xrange(self._dim):
					if self._partial_var_active_set[d] <= self._tol_dims[d]:

						self._levels_converged_dims[d] = curr_max_levels[d]
						self._levels_converged_dims_all.append(self._levels_converged_dims.copy())
			else:
				global_min_levels = np.amin(self.levels_converged_dims_all, axis=0)

				for d in xrange(self._dim):
					if self._partial_var_active_set[d] <= self._tol_dims[d]:
						if curr_max_levels[d] < global_min_levels[d]:
							self._levels_converged_dims[d] = curr_max_levels[d]
							self._levels_converged_dims_all.append(self._levels_converged_dims.copy())
						else:
							self._levels_converged_dims[d] = global_min_levels[d]

	def check_termination_criterion(self):

		max_level = np.max(self._multiindex_set)
		if len(self._A.values()) == 0 or self._eta <= self._tol or max_level >= self._max_level or \
									self.__is_less_or_equal_then(self._partial_var_active_set, self._tol_dims):
			self._stop_adaption = True 

	def serialize_data(self, serialization_file):
		
		with open(serialization_file, "wb") as output_file:
			data = [self._key_O, self._O, self._key_A, self._A, self._key_local_error, self._local_error, \
					self._eta, self._multiindex_set, self._levels_converged_dims, self._levels_converged_dims_all, \
					self._partial_var_active_set, self._local_basis_local, self._local_basis_global]
			dump(data, output_file)

		output_file.close()

	def unserialize_data(self, serialization_file):

		with open(serialization_file, "rb") as input_file:
			self._key_O, self._O, self._key_A, self._A, self._key_local_error, self._local_error, \
					self._eta, self._multiindex_set, self._levels_converged_dims, self._levels_converged_dims_all, \
					self._partial_var_active_set, self._local_basis_local, self._local_basis_global = load(input_file)

		input_file.close()
