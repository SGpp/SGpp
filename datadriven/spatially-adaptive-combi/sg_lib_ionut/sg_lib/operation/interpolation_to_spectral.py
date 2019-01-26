from onedim import *
from abstract_operation import *
from ast import literal_eval

class InterpolationToSpectral(AbstractOperation):
	def __init__(self, dim, linear_growth_factor, left_bounds, right_bounds, weights):
		
		self._dim 					= dim
		self.__linear_growth_factor = linear_growth_factor
		self.__left_bounds 			= left_bounds
		self.__right_bounds 		= right_bounds
		self.__weights 				= weights
		
		self._sg_func_evals_key	= 0
		self._sg_func_evals_lut = OrderedDict()
		self._fg_func_evals 	= None

		self._grid_points_1D       	= []
		self._global_indices_dict   = OrderedDict()
		self._no_fg_grid_points     = 0

		self.__local_basis 	= None
		self.__global_basis = None

	def __eval_ND_orth_poly(self, degs, x):

		poly_eval = 1.
		for d in xrange(self._dim):
			poly_eval *= eval_1D_orth_poly(degs[d], self.__left_bounds[d], self.__right_bounds[d], self.__weights[d], x[d])

		return poly_eval

	def __get_no_1D_grid_points(self, level):

		no_points = 0

		if self.__linear_growth_factor == 'sym':
			no_points = 2*level - 1
		else:
			if self.__linear_growth_factor == 1:
				no_points = level
			elif self.__linear_growth_factor >= 2:
				if level == 1:
					no_points = 1
				else:
					no_points = self.__linear_growth_factor*level - 1
			else:
				raise NotImplementedError

		return no_points

	def __get_orth_poly_degs_inf_norm(self, multiindex):

		if np.sum(multiindex) == self._dim:
			tensorized_degrees = [[0 for i in xrange(self._dim)]]
		else:
			tensorized_degrees = self.__local_basis[repr(multiindex.tolist())]

		return tensorized_degrees

	def __get_orth_poly_basis_global(self, multiindex_set):

	 	orth_poly_basis_global = self.__global_basis

	 	return orth_poly_basis_global

	def __get_orth_poly_basis_active_set(self, active_set):

		active_set_basis = []

		for multiindex in active_set:
			degs = self.__get_orth_poly_degs_inf_norm(multiindex)
			
			for deg in degs:
				if deg not in active_set_basis:
					active_set_basis.append(deg)

		active_set_basis = np.array(active_set_basis, dtype=int)

		return active_set_basis

	def __get_spectral_coeff_local_dict(self, func_evals, multiindex, tensorized_degrees):

		spectral_coeff_fg = self.__get_spectral_coeff_local(func_evals, multiindex)

		spectral_coeff_dict = OrderedDict()
		for degrees, coeff_fg in zip(tensorized_degrees, spectral_coeff_fg):
			spectral_coeff_dict[repr(degrees)] = coeff_fg

		return spectral_coeff_dict

	def __get_spectral_coeff_global(self, func_evals, multiindex, orth_poly_basis):

		spectral_coeff 			= np.zeros(len(orth_poly_basis))
		tensorized_degrees 		= self.__get_orth_poly_degs_inf_norm(multiindex)
		curr_spectral_coeff 	= self.__get_spectral_coeff_local_dict(func_evals, multiindex, tensorized_degrees)

		if type(orth_poly_basis) is np.ndarray:
			orth_poly_basis = orth_poly_basis.tolist()
		
		for local_pos in curr_spectral_coeff.keys():
			local_pos_list = list(literal_eval(local_pos))

			if local_pos_list in orth_poly_basis:
				index 					= orth_poly_basis.index(local_pos_list)
				spectral_coeff[index] 	= curr_spectral_coeff[local_pos]

		return spectral_coeff

	def __get_spectral_coeff_delta_dict(self, curr_multiindex, multiindex_set):
		
		spectral_coeff_delta 	= self.get_spectral_coeff_delta(curr_multiindex, multiindex_set)
		orth_poly_basis_local 	= self.__get_orth_poly_degs_inf_norm(curr_multiindex)

		spectral_coeff_dict = OrderedDict()
		for degrees, coeff_fg in zip(orth_poly_basis_local, spectral_coeff_delta):
			spectral_coeff_dict[repr(degrees)] = coeff_fg

		return spectral_coeff_dict

	def __get_spectral_coeff_active_set(self, curr_multiindex, multiindex_set, orth_poly_basis):

		spectral_coeff 			= np.zeros(len(orth_poly_basis))
		curr_spectral_coeff 	= self.__get_spectral_coeff_delta_dict(curr_multiindex, multiindex_set)
		orth_poly_basis 		= orth_poly_basis.tolist()

		for local_pos in curr_spectral_coeff.keys():
			local_pos_list = list(literal_eval(local_pos))

			if local_pos_list in orth_poly_basis:
				index 					= orth_poly_basis.index(local_pos_list)
				spectral_coeff[index] 	= curr_spectral_coeff[local_pos]

		return spectral_coeff

	def _eval_operation_fg(self, curr_func_evals, multiindex, x):

		interp_fg 	 	= 0.
		poly_eval_all 	= []

		for d in xrange(self._dim):
			index 				= multiindex[d] 
			grid_1D 			= self._grid_points_1D[d][index - 1]
			barycentric_weights	= get_1D_barycentric_weights(grid_1D)

			poly_eval = []
			for j in xrange(len(grid_1D)):
				p_eval = eval_1D_barycentric_interpolant(grid_1D, barycentric_weights, j, x[d])
				poly_eval.append(p_eval)

			poly_eval_all.append(poly_eval)

		tensorized_basis_val = np.array([np.prod(interp_pair) for interp_pair in list(product(*poly_eval_all))], dtype=np.float64)

		for i in xrange(len(tensorized_basis_val)):
			interp_fg += tensorized_basis_val[i]*curr_func_evals[i]

		return interp_fg

	# map interpolation to spectral projection
	def __get_spectral_coeff_local(self, curr_func_evals, multiindex):

		orth_multiindex_degs = self.__get_orth_poly_degs_inf_norm(multiindex)

		curr_1D_points = []
		for d in xrange(self._dim):
			index 	= multiindex[d] 
			grid_1D = self._grid_points_1D[d][index - 1]
			curr_1D_points.append(grid_1D)

		tensorized_leja_points = np.array([leja_points for leja_points in list(product(*curr_1D_points))], dtype=np.float64)

		assert len(orth_multiindex_degs) == len(tensorized_leja_points)

		rhs = np.zeros(len(orth_multiindex_degs))
		for i, leja_point in enumerate(tensorized_leja_points):
			rhs[i] = self._eval_operation_fg(curr_func_evals, multiindex, leja_point)

		orth_basis_matrix = np.zeros((len(orth_multiindex_degs), len(orth_multiindex_degs)))
		for i, leja_point in enumerate(tensorized_leja_points):
			for j, poly_degs in enumerate(orth_multiindex_degs):
				orth_basis_matrix[i, j] = self.__eval_ND_orth_poly(poly_degs, leja_point)

		spectral_coeff_fg = np.linalg.solve(orth_basis_matrix, rhs)

		return spectral_coeff_fg

	def get_local_global_basis(self, adaptivity_obj):

		self.__local_basis 	= adaptivity_obj.local_basis_local
		self.__global_basis = adaptivity_obj.local_basis_global

	def get_spectral_coeff_delta(self, curr_multiindex, multiindex_set):
		
		multiindex_dict 		= self._get_multiindex_dict(multiindex_set)
		orth_poly_basis_local 	= np.array(self.__get_orth_poly_degs_inf_norm(curr_multiindex), dtype=int)
		no_spectral_coeff 		= len(orth_poly_basis_local)
		spectral_coeff_delta 	= np.zeros(no_spectral_coeff)

		differences_indices, differences_signs 	= self._get_differences_sign(curr_multiindex)

		keys_differences 	= differences_indices.keys()
		keys_signs 			= differences_signs.keys() 
		
		for key in keys_differences:
			differences 	= differences_indices[key]
			curr_func_evals = self._fg_func_evals[multiindex_dict[repr(differences.tolist())]]
			
			sign 					= differences_signs[key] 
			spectral_coeff_delta 	+= self.__get_spectral_coeff_global(curr_func_evals, \
																	differences, orth_poly_basis_local)*sign

		return spectral_coeff_delta			

	def get_spectral_coeff_sg(self, multiindex_set):

		multiindex_dict 		= self._get_multiindex_dict(multiindex_set)
		orth_poly_basis_global 	= self.__get_orth_poly_basis_global(multiindex_set)
		no_spectral_coeff 		= len(orth_poly_basis_global)
		spectral_coeff 			= np.zeros(no_spectral_coeff)

		for index, multiindex in enumerate(multiindex_set):
			differences_indices, differences_signs = self._get_differences_sign(multiindex)

			keys_differences 	= differences_indices.keys()
			keys_signs 			= differences_signs.keys() 

			for key in keys_differences:
				differences 	= differences_indices[key]
				curr_func_evals = self._fg_func_evals[multiindex_dict[repr(differences.tolist())]]
				sign 			= differences_signs[key] 

				curr_spectral_coeff = self.__get_spectral_coeff_global(curr_func_evals, differences, orth_poly_basis_global)
				spectral_coeff 		+= sign*curr_spectral_coeff

		return spectral_coeff

	def eval_operation_delta(self, curr_multiindex, multiindex_set, x):

		interp_delta = 0.
		
		multiindex_dict 						= self._get_multiindex_dict(multiindex_set)
		differences_indices, differences_signs 	= self._get_differences_sign(curr_multiindex)

		keys_differences 	= differences_indices.keys()
		keys_signs 			= differences_signs.keys() 

		for key in keys_differences:
			differences 	= differences_indices[key]
			curr_func_evals = self._fg_func_evals[multiindex_dict[repr(differences.tolist())]]

			sign 			= differences_signs[key] 
			interp_delta 	+= self._eval_operation_fg(curr_func_evals, differences, x)*sign

		return interp_delta

	def eval_operation_sg(self, multiindex_set, x):

		interp_sg = 0.

		for multiindex in multiindex_set:
			interp_delta = self.eval_operation_delta(multiindex, multiindex_set, x) 
			interp_sg 	+= interp_delta

		return interp_sg

	def get_mean(self, spectral_coeff):
		
		mean = spectral_coeff[0]

		return mean

	def get_variance(self, spectral_coeff):
		
		var = np.sum([spectral_coeff[i]**2 for i in xrange(1, len(spectral_coeff))])

		return var

	def get_directional_var_active_set(self, global_multiindex_set, active_set):
		
		directional_var = np.zeros(self._dim)

		orth_poly_basis_active 	= self.__get_orth_poly_basis_active_set(active_set)
		spectral_coeff 			= np.zeros(len(orth_poly_basis_active))

		for multiindex in active_set:
			spectral_coeff_active_set 	= self.__get_spectral_coeff_active_set(multiindex, global_multiindex_set, orth_poly_basis_active)
			spectral_coeff 				+= spectral_coeff_active_set

		for d in xrange(self._dim):
			
			mindex = []
			for i in xrange(len(orth_poly_basis_active)):
				if orth_poly_basis_active[i][d] != 0:
					mindex.append(i)

 			directional_var[d] = np.sum([spectral_coeff[j]**2 for j in mindex])

		return directional_var

	def get_local_var_level_active_set(self, global_multiindex_set, multiindex):
		
		directional_var = np.zeros(self._dim + 1)

		orth_poly_basis_local 	= self.__get_orth_poly_degs_inf_norm(multiindex)
		spectral_coeff_local 	= self.get_spectral_coeff_delta(multiindex, global_multiindex_set)

		mindex_interact = []
		for d in xrange(self._dim):
			
			mindex_main = []
			for i in xrange(len(orth_poly_basis_local)):
				if orth_poly_basis_local[i][d] != 0 and np.count_nonzero(orth_poly_basis_local[i]) == 1:
					mindex_main.append(i)
				elif orth_poly_basis_local[i][d] != 0 and np.count_nonzero(orth_poly_basis_local[i]) >= 2:
					mindex_interact.append(i)

 			directional_var[d] = np.sum([spectral_coeff_local[j]**2 for j in mindex_main])

 		directional_var[self._dim] = np.sum([spectral_coeff_local[j]**2 for j in mindex_interact])

		return directional_var

	def get_total_var_level_active_set(self, global_multiindex_set, multiindex):
		
		spectral_coeff_local = self.get_spectral_coeff_delta(multiindex, global_multiindex_set)

		total_var_level = np.sum(spectral_coeff_local[1:]**2)

		return total_var_level

	def get_first_order_sobol_indices(self, spectral_coeff, multiindex_set):
		
		sobol_indices = np.zeros(self._dim)

		orth_poly_basis_global 	= self.__get_orth_poly_basis_global(multiindex_set)
		Var 					= np.sum([spectral_coeff[i]**2 for i in xrange(1, len(spectral_coeff))])

		for d in xrange(self._dim):
			
			mindex = []
			for i in xrange(len(orth_poly_basis_global)):
				if orth_poly_basis_global[i][d] != 0 and np.count_nonzero(orth_poly_basis_global[i]) == 1:
					mindex.append(i)

 			sobol_indices[d] = np.sum([spectral_coeff[j]**2 for j in mindex])/Var

		return sobol_indices

	def get_total_sobol_indices(self, spectral_coeff, multiindex_set):
		
		sobol_indices = np.zeros(self._dim)

		orth_poly_basis_global 	= self.__get_orth_poly_basis_global(multiindex_set)
		Var 					= np.sum([spectral_coeff[i]**2 for i in xrange(1, len(spectral_coeff))])

		for d in xrange(self._dim):
			
			mindex = []
			for i in xrange(len(orth_poly_basis_global)):
				if orth_poly_basis_global[i][d] != 0:
					mindex.append(i)

 			sobol_indices[d] = np.sum([spectral_coeff[j]**2 for j in mindex])/Var

		return sobol_indices

	def serialize_results(self, E, Var, Sobol_indices, serialization_file):
	    
	    with open(serialization_file, "wb") as output_file:
	    	data = [E, Var, Sobol_indices]
	        dump(data, output_file)

	    output_file.close()

	def unserialize_results(self, serialization_file):

	    with open(serialization_file, "rb") as input_file:
	    	E, Var, Sobol_indices = load(input_file)
	           
	    input_file.close() 

	    return E, Var, Sobol_indices
