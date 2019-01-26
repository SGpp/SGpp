from abstract_adapt_operation import *

class SpectralScores(DimensionAdaptivity):

	def __init__(self, dim, tols, init_multiindex, max_level, level_to_nodes, spectral_op_obj):
		
		self._dim 				= dim
		self._tols 				= tols
		self._init_multiindex 	= init_multiindex
		self._max_level 		= max_level
		self._level_to_nodes 	= level_to_nodes

		self.__spectral_op_obj = spectral_op_obj

		self._eta = 0.

		self._O 				= OrderedDict()
		self._A 				= OrderedDict() 
		self._local_error 		= OrderedDict()
		self._key_O 			= 0
		self._key_A 			= 0
		self._key_local_error 	= 0

		self._multiindex_set 	= []
		self._init_no_points 	= 0

		self._stop_adaption = False

		self._local_basis_global 	= None
		self._local_basis_local 	= OrderedDict()

	def __get_local_score(self, curr_multiindex):
		
		local_score 	= 0.
		local_variances = self.__spectral_op_obj.get_local_var_level_active_set(self._multiindex_set, curr_multiindex)

		if np.sum(curr_multiindex) == self._dim:
			local_score = 1
		else:
			for d in xrange(self._dim + 1):
				if local_variances[d] >= self._tols[d]:
					local_score += 1

		return local_score

	def __get_max_score(self):

		max_scores_pos 	= np.where(np.array([local_score == np.amax(self._local_error.values()) for local_score in self._local_error.values()]))[0]
		max_scores_keys = np.array([self._local_error.keys()[max_pos] for max_pos in max_scores_pos])

		max_elem 	= 0.
		max_key 	= 0
		if len(max_scores_keys) >= 2:
			for key in max_scores_keys:
				multiindex 		= self._A[key]

				local_variance 	= self.__spectral_op_obj.get_total_var_level_active_set(self._multiindex_set, multiindex)

				if local_variance >= max_elem:
					max_elem 	= local_variance
					max_key 	= key

		else:
			max_key = max_scores_keys[0]

		return max_key

	def init_adaption(self):

		self._multiindex_set.append(self._init_multiindex)

		self._key_O 								= -1
		self._key_A 								= 0
		self._key_local_error 						= 0
		self._A[self._key_A] 						= self._init_multiindex
		local_score 								= self.__get_local_score(self._init_multiindex)
		self._local_error[self._key_local_error] 	= local_score
		self._eta 	 								= local_score

		self._local_basis_local[repr(self._init_multiindex)] 	= self._get_local_hierarchical_basis(self._init_multiindex)
		self._local_basis_global 								= self._get_local_hierarchical_basis(self._init_multiindex)

	def do_one_adaption_step_preproc(self):

		local_multiindices = []

		max_index 	= self.__get_max_score()
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

				local_basis_neighbor = np.array([self._get_no_1D_grid_points(n) - 1 for n in neighbor], dtype=int)
				self._update_local_basis(neighbor.tolist(), local_basis_neighbor)

		local_multiindices = np.array(local_multiindices, dtype=int)

		return local_multiindices

	def do_one_adaption_step_postproc(self, curr_multiindices):

		for multiindex in curr_multiindices:
			self._key_local_error 					+= 1
			local_score					 			 = self.__get_local_score(multiindex)
			self._local_error[self._key_local_error] = local_score

			self._eta += local_score

	def check_termination_criterion(self):

		max_level = np.max(self._multiindex_set)
		if len(self._A.values()) == 0 or np.sum(self._local_error.values()) == 0. or max_level >= self._max_level:
			self._stop_adaption = True

	def serialize_data(self, serialization_file):
		
		with open(serialization_file, "wb") as output_file:
			data = [self._key_O, self._O, self._key_A, self._A, self._key_local_error, self._local_error, self._multiindex_set, \
							self._local_basis_local, self._local_basis_global]
			dump(data, output_file)

		output_file.close()

	def unserialize_data(self, serialization_file):

		with open(serialization_file, "rb") as input_file:
			self._key_O, self._O, self._key_A, self._A, self._key_local_error, self._local_error, self._multiindex_set, \
							self._local_basis_local, self._local_basis_global = load(input_file)

		input_file.close()