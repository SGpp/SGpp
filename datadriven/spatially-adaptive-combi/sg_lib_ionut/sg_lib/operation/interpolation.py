from onedim import *
from abstract_operation import *

class Interpolation(AbstractOperation):
	def __init__(self, dim):
		
		self._dim = dim

		self._sg_func_evals_key	= 0
		self._sg_func_evals_lut = OrderedDict()
		self._fg_func_evals 	= None

		self._grid_points_1D       	= []
		self._global_indices_dict   = OrderedDict()
		self._no_fg_grid_points     = 0

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