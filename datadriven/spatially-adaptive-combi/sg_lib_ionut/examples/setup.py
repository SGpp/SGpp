import numpy as np

dim 			= 3
left_bounds 	= np.zeros(dim)
right_bounds 	= np.ones(dim)
grid_level 		= 1
growth_factor 	= 2
ref_step 		= 0
total_no_gp 	= 0

tol 		= 1e-4
tol_dims 	= 1e-5*np.ones(dim)
max_level 	= 10

weight = lambda x: 1.0

weights 	= []
for d in xrange(dim):
	weights.append(weight)

init_multiindex = np.ones(dim, dtype=int)

eval_point = 0.33*np.ones(dim)

data_file 	= 'serialized_files/data.pickle'
ref_file 	= 'serialized_files/refinement.pickle'
info_file 	= 'serialized_files/ref_info.pickle'
mindex_file = 'serialized_files/new_multiindices.pickle'