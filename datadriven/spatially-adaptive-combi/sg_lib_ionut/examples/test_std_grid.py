from sg_op.algebraic.multiindex import *
from matplotlib.pyplot import *
from matplotlib.patches import *

if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 7
	growth_factor 	= 2

	weight_dim1 = lambda x: 1.0
	weight_dim2 = lambda x: 1.0

	weights 	= []
	weights.append(weight_dim1)
	weights.append(weight_dim2)

	Multiindex_obj 	= Multiindex(dim)

	multiindex_set 	= Multiindex_obj.get_std_total_degree_mindex(grid_level)

	# plotting
	fig = figure()
	ax 	= fig.add_subplot(111, aspect='equal')

	temp 				= np.ones(dim, dtype=int)
	max_reached_level 	= grid_level

	for p in [
		Rectangle(
	    tuple(multiindex - temp),
	    1.0,
	    1.0,
	    hatch='\\',
	    facecolor='lightblue',
	    edgecolor='black', 
	    linewidth=4.0,
	    label='standard multiindex set' ) for multiindex in multiindex_set]:
		ax.add_patch(p)
		old_index_set_handle = p

	for multiindex in multiindex_set:
		x_coord = (2*tuple(multiindex - temp)[0] + 0.5)/2.
		y_coord = (2*tuple(multiindex - temp)[1] + 0.9)/2.

		ax.text(x_coord, y_coord, str(tuple(multiindex)), fontsize=20)

	
	legend(handles=[old_index_set_handle], loc='best', fontsize=25)

	x_coord = (2*max_reached_level - 3)/2.
	y_coord = (2*max_reached_level - 2)/2.	

	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_xlim([0, max_reached_level])
	ax.set_ylim([0, max_reached_level])

	show()