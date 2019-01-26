from onedim import *
from abstract_operation import *
from ast import literal_eval

class SpectralProjection(AbstractOperation):
	def __init__(self, dim, left_bounds, right_bounds, weights):
		
		self._dim 			= dim
		self.left_bounds 	= left_bounds
		self.right_bounds 	= right_bounds
		self.weights 		= weights

		self._sg_func_evals_key	= 0
		self._sg_func_evals_lut = OrderedDict()
		self._fg_func_evals 	= None

		self._grid_points_1D       	= []
		self._global_indices_dict   = OrderedDict()
		self._no_fg_grid_points     = 0

	def __get_orth_poly_degs_fg(self, multiindex):

		degrees_all = []

		for d in xrange(self._dim):
			grid_1D_len = 2*multiindex[d] - 1
			pl 			= (grid_1D_len - 1)/2
			
			degrees__dim_d = []
			for p in xrange(pl + 1):
				degrees__dim_d.append(p)
				
			degrees_all.append(degrees__dim_d)

		tensorized_degrees = list(product(*degrees_all))

		return tensorized_degrees

	# def __get_orth_poly_basis_global(self, multiindex_set):

	# 	max_level_deg 			= (2*np.max(multiindex_set) - 2)/2
	# 	orth_poly_basis_global 	= Multiindex(self._dim).get_poly_mindex(max_level_deg)

	# 	return orth_poly_basis_global

	def __get_orth_poly_basis_global(self, multiindex_set):

		if type(multiindex_set) is not np.ndarray:
			multiindex_set = np.array(multiindex_set, dtype=int)

		max_degs = np.zeros(self._dim, dtype=int)

		for d in xrange(self._dim):
			max_degs[d] = (2*np.max(multiindex_set.T[d]) - 2)/2

		max_level_deg 		= (2*np.max(multiindex_set) - 2)/2
		tensorized_degrees 	= Multiindex(self._dim).get_poly_mindex(max_level_deg).tolist()

		g_basis = []
	
		for deg in tensorized_degrees:
			if deg not in g_basis and np.sum(deg) <= np.max(max_degs):
				
				is_feasible = 1.
				for d in xrange(self._dim):
					if deg[d] > max_degs[d]:
						is_feasible = 0.
						break

				if is_feasible:
					g_basis.append(deg)

		g_basis = np.array(g_basis, dtype=int)

		return g_basis

	def __get_orth_poly_basis_active_set(self, active_set):

		active_set_basis = []

		for multiindex in active_set:
			degs = self.__get_orth_poly_degs_fg(multiindex)
			
			for deg in degs:
				if deg not in active_set_basis:
					active_set_basis.append(deg)

		active_set_basis = np.array(active_set_basis, dtype=int)

		return active_set_basis

	def __get_spectral_coeff_local(self, func_evals, multiindex):

		orth_poly_all 		= []
		quad_weights_all 	= []

		for d in xrange(self._dim):
			index 				= multiindex[d] 
			grid_1D 			= self._grid_points_1D[d][index - 1]
			barycentric_weights	= get_1D_barycentric_weights(grid_1D)

			quad_weights = compute_1D_quad_weights(grid_1D, self.left_bounds[d], self.right_bounds[d], self.weights[d])
			quad_weights_all.append(quad_weights)

			pl = (len(grid_1D) - 1)/2
			
			orth_poly__dim_d_all = []
			for p in xrange(pl + 1):

				orth_poly__dim_d = []
				for q in xrange(len(grid_1D)):
					orth_poly = eval_1D_orth_poly(p, self.left_bounds[d], self.right_bounds[d], self.weights[d], grid_1D[q])
					orth_poly__dim_d.append(orth_poly)

				orth_poly__dim_d_all.append(orth_poly__dim_d)			
			orth_poly_all.append(orth_poly__dim_d_all)

		tensorized_weights 		= [np.prod(weights_pair) for weights_pair in list(product(*quad_weights_all))]
		tensorized_poly_bases 	= [[np.prod(poly_pair) for poly_pair in list(product(*poly_pairs))] for poly_pairs in list(product(*orth_poly_all))]

		spectral_coeff_fg = np.zeros(len(tensorized_poly_bases))
		for p, tensorized_poly_basis in enumerate(tensorized_poly_bases):
			for q in xrange(len(tensorized_poly_basis)):
				spectral_coeff_fg[p] += func_evals[q]*tensorized_poly_basis[q]*tensorized_weights[q]

		return spectral_coeff_fg

	def __get_spectral_coeff_local_dict(self, func_evals, multiindex, tensorized_degrees):

		spectral_coeff_fg = self.__get_spectral_coeff_local(func_evals, multiindex)

		spectral_coeff_dict = OrderedDict()
		for degrees, coeff_fg in zip(tensorized_degrees, spectral_coeff_fg):
			spectral_coeff_dict[repr(degrees)] = coeff_fg

		return spectral_coeff_dict

	def __get_spectral_coeff_global(self, func_evals, multiindex, orth_poly_basis):

		spectral_coeff 			= np.zeros(len(orth_poly_basis))
		tensorized_degrees 		= self.__get_orth_poly_degs_fg(multiindex)
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
		orth_poly_basis_local 	= self.__get_orth_poly_degs_fg(curr_multiindex)

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
		
		spectral_fg = 0.

		tensorized_degrees 		= self.__get_orth_poly_degs_fg(multiindex)
		spectral_coeff_fg 		= self.__get_spectral_coeff_local(curr_func_evals, multiindex)
		poly_eval_all 			= []

		for p in tensorized_degrees:

			poly_eval__dim_d = []
			for d in xrange(self._dim):
				poly_eval = eval_1D_orth_poly(p[d], self.left_bounds[d], self.right_bounds[d], self.weights[d], x[d])
				poly_eval__dim_d.append(poly_eval)

			poly_eval_all.append(poly_eval__dim_d)

		tensorized_poly_eval = np.array([np.prod(poly_pair) for poly_pair in poly_eval_all])

		for i, basis_val in enumerate(tensorized_poly_eval):
			spectral_fg += spectral_coeff_fg[i] * basis_val

		return spectral_fg

	def get_spectral_coeff_delta(self, curr_multiindex, multiindex_set):
		
		multiindex_dict 		= self._get_multiindex_dict(multiindex_set)
		orth_poly_basis_local 	= np.array(self.__get_orth_poly_degs_fg(curr_multiindex), dtype=int)
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
		
		spectral_delta = 0.

		multiindex_dict 						= self._get_multiindex_dict(multiindex_set)
		differences_indices, differences_signs 	= self._get_differences_sign(curr_multiindex)

		keys_differences 	= differences_indices.keys()
		keys_signs 			= differences_signs.keys() 

		for key in keys_differences:
			differences 	= differences_indices[key]
			curr_func_evals = self._fg_func_evals[multiindex_dict[repr(differences.tolist())]]
			
			sign 			= differences_signs[key] 
			spectral_delta 	+= self._eval_operation_fg(curr_func_evals, differences, x)*sign

		return spectral_delta

	def eval_operation_sg(self, multiindex_set, x):

		spectral_sg = 0.

		for multiindex in multiindex_set:
			spectral_delta 	= self.eval_operation_delta(multiindex, multiindex_set, x) 
			spectral_sg 	+= spectral_delta

		return spectral_sg

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
