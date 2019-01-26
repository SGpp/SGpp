from onedim import *
from abstract_operation import *

class Quadrature(AbstractOperation):
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

	def _eval_operation_fg(self, curr_func_evals, multiindex):

		quad_fg 	 		= 0.
		quad_weights_all 	= []

		for d in xrange(self._dim):
			index 				= multiindex[d] 
			grid_1D 			= self._grid_points_1D[d][index - 1]

			quad_weights = compute_1D_quad_weights(grid_1D, self.left_bounds[d], self.right_bounds[d], self.weights[d])
			quad_weights_all.append(quad_weights)

		tensorized_weights = np.array([np.prod(weights_pair) for weights_pair in list(product(*quad_weights_all))], dtype=np.float64)

		for i in xrange(len(tensorized_weights)):
			quad_fg += tensorized_weights[i]*curr_func_evals[i]

		return quad_fg

	def eval_operation_delta(self, curr_multiindex, multiindex_set):

		quad_delta = 0.
		
		multiindex_dict 						= self._get_multiindex_dict(multiindex_set)
		differences_indices, differences_signs 	= self._get_differences_sign(curr_multiindex)

		keys_differences 	= differences_indices.keys()
		keys_signs 			= differences_signs.keys() 

		for key in keys_differences:
			differences 	= differences_indices[key]
			curr_func_evals = self._fg_func_evals[multiindex_dict[repr(differences.tolist())]]

			sign 		= differences_signs[key] 
			quad_delta 	+= self._eval_operation_fg(curr_func_evals, differences)*sign

		return quad_delta	

	def eval_operation_sg(self, multiindex_set):

		quad_sg = 0.

		for multiindex in multiindex_set:
			quad_delta 	= self.eval_operation_delta(multiindex, multiindex_set) 
			quad_sg 	+= quad_delta

		return quad_sg