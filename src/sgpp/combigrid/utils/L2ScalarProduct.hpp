/*
 * L2ScalarProduct.hpp
 *
 *  Created on: 25.02.2012
 *      Author: ckow
 */

#ifndef L2SCALARPRODUCT_HPP_
#define L2SCALARPRODUCT_HPP_
#include <math.h>
#include <vector>

namespace combigrid {

using namespace std;
template<typename ELEMENT>
class L2ScalarProduct {

public:

	L2ScalarProduct(int dim, vector<bool> hasBoundary, vector<int> gridShape) {
		dim_ = dim;
//		boundaryPoints_=hasBoundary;
//		gridShape_=gridShape;
		boundaryPoints_ = new bool[dim];
		levels_ = new int[dim];
		numPts_ = 1;
		for (int i = 0; i < dim; ++i) {
			boundaryPoints_[i] = hasBoundary[i];
			levels_[i] = gridShape[i];
			if (boundaryPoints_[i])
				numPts_ *= pow(2, levels_[i]) + 1;
			else
				numPts_ *= pow(2, levels_[i]) + 1;

		}

		//calculate gram values
		gramValues = calculateGramValues();
	}

	//template<typename ELEMENT>
	//L2ScalarProduct<ELEMENT>::~L2ScalarProduct() {
	//	// TODO Auto-generated destructor stub
	//}

	ELEMENT return_l2_scalar_product(vector<ELEMENT> u_in,
			vector<ELEMENT> v_in) {
		ELEMENT* u = new ELEMENT[u_in.size()];
		ELEMENT* v = new ELEMENT[v_in.size()];

		for (int i = 0; i < u_in.size(); ++i) {
			u[i] = u_in[i];
			v[i] = v_in[i];
//			cout<<u[i]<<'\t'<<v[i]<<endl;
		}

		int dimension = dim_;
		int* levels = levels_;
		int* grid_shape = new int[dimension];
		int* grid_strides = new int[dimension];
		int* gridpt_shape = new int[dimension];
		int* gridpt_strides = new int[dimension];

		int num_grids;
		int num_gridpts;

		// Get the shape and strides of the grid and grid points given levels an dimension
		get_grid_info(levels, dimension, grid_shape, grid_strides, num_grids,
				gridpt_shape, gridpt_strides, num_gridpts);
		// Get the corner shifts needed to retrieve the corners from each grid element
		int num_corners = pow(2, dimension);
		int* corner_shifts = new int[num_corners];
		get_corner_shifts(gridpt_shape, gridpt_strides, dimension,
				corner_shifts);

		// Get the gram values of the basis for the grid (that is the inner product of each basis with every other basis)
		//	ELEMENT * all_gram_values = calculateGramValues();
		double * all_gram_values = new double[num_corners*num_corners];
		for (int i = 0; i < num_corners*num_corners; ++i) {
			all_gram_values[i]=gramValues[i];
		}

		ELEMENT* temp_u = new ELEMENT[num_corners];
		ELEMENT* temp_v = new ELEMENT[num_corners];

		double* agv_p = all_gram_values;

		// Uses the gram values to calculate the scalar product
		ELEMENT S = 0;

		int dim_counter = dimension;

		_return_l2_scalar_product_(dim_counter, grid_shape, gridpt_strides, u,
				v, all_gram_values, corner_shifts, num_corners, temp_u, temp_v,
				S);

		// delete all temporaries
		delete gridpt_shape;
		delete gridpt_strides;
		delete grid_shape;
		delete grid_strides;

		delete temp_u;
		delete temp_v;

		delete all_gram_values;
		delete corner_shifts;

		return S;
	}
	/**
	 Return the all gram values (that is, the scalar product of each basis with every other basis)

	 The output format is described as follows:
	 If in the case of of one dimension, there are basically two basis phi_0, phi_1,
	 the output is then an array of phi_0*phi_0, phi_0*phi_1, phi_1*phi_0, phi_1*phi_1.

	 For efficiency, phi_1*phi_0, phi_1*phi_1 is not actually calculated again, but permutated from
	 phi_0*phi_0, phi_0*phi_1, that is why reflect_vector is used.

	 **/

	double * calculateGramValues() {
		int dimension = dim_;
		int num_corners = int(pow(2., dimension));
		int n = int(pow(2., levels_[0]));
		double* gram_at_level = new double[2];
		gram_at_level[0] = 1 / (3.0 * n);
		gram_at_level[1] = 1 / (6.0 * n);

		double* gram_values = new double[num_corners];

		copy(gram_values, gram_at_level, 2);

		if (dimension > 1) {

			int gram_size = 2;

			double* gram_temp = new double[num_corners];

			for (int i = 1; i < dimension; i++) {

				n = int(pow(2., levels_[i]));

				gram_at_level[0] = 1 / (3.0 * n);
				gram_at_level[1] = 1 / (6.0 * n);

				getOuterProduct(gram_at_level, 2, gram_values, gram_size,
						gram_temp);

				gram_size *= 2;
				copy(gram_values, gram_temp, gram_size);

			}

			delete gram_temp;

		}

		delete gram_at_level;

		//Reflects
		int* choice = new int[dimension];
		setZero(choice, dimension);

		int* choice_shape = new int[dimension];
		setConstant(choice_shape, dimension, 2);

		int* gram_shape = new int[dimension];
		setConstant(gram_shape, dimension, 2);

		int* gram_strides = new int[dimension];
		gram_strides[0] = 1;
		for (int i = 1; i < dimension; i++)
			gram_strides[i] = gram_strides[i - 1] * gram_shape[i - 1];

		double* all_gram_values = new double[num_corners * num_corners];

		double* agv_p = all_gram_values;

		for (int i = 0; i < num_corners; i++, agv_p += num_corners) {
			getReflectVector(gram_strides, gram_shape, dimension, choice,
					gram_values, agv_p);

			getNextCoordinate(choice_shape, choice, dimension);

		}
		agv_p = all_gram_values;

		delete choice;
		delete choice_shape;
		delete gram_shape;
		delete gram_strides;
		delete gram_values;

		return all_gram_values;
	}

	void getOuterProduct(double*& vec1, int vec1_size, double*& vec2,
			int vec2_size, double*& output) {
		double* output_p = output;
		double* vec1_p = vec1;
		double* vec2_p = vec2;

		for (int i = 0; i < vec1_size; i++, vec1_p++) {
			for (int j = 0; j < vec2_size; j++, vec2_p++) {
				*output_p = *vec1_p * *vec2_p;

				output_p++;
			}

			vec2_p = vec2;
		}

	}

	void setZero(int * vec, int size) {
		int* vec_pointer = vec;
		for (int i = 0; i < size; i++, vec_pointer++)
			*vec_pointer = 0.;
	}

	void setConstant(int *vec, int size, int constant) {
		int* vec_pointer = vec;
		for (int i = 0; i < size; i++, vec_pointer++)
			*vec_pointer = constant;
	}

	//int* double*
	void _getReflectVector_(int*& strides, int*& shape, int& dimension,
			int *&axis, double*& vec_pointer, double*& reflected_vec_pointer) { /*
			 Reflects the array given by pointer vec with strides and shape against axis and write on to reflected_vec.

			 axis must be of the form (0,1,0...,1) etc where 0 means no reflection and 1 means reflection
			 a temporary int initialized to dimension must be supplied.
			 */
		if (dimension == 1) {
			if (axis[0] == 0) {
				for (int i = 0; i < shape[0]; i++) {

					*reflected_vec_pointer = *vec_pointer;
					vec_pointer++;
					reflected_vec_pointer++;

				}
				vec_pointer -= strides[1];
			} else {
				for (int i = 0; i < shape[0]; i++) {

					*reflected_vec_pointer = *vec_pointer;
					vec_pointer--;
					reflected_vec_pointer++;

				}
				vec_pointer += strides[1];
			}

		} else {
			dimension--;
			if (axis[dimension] == 0) {
				for (int i = 0; i < shape[dimension]; i++) {
					_getReflectVector_(strides, shape, dimension, axis,
							vec_pointer, reflected_vec_pointer);
					vec_pointer += strides[dimension];
				}

				vec_pointer -= strides[dimension + 1];
			} else {

				for (int i = 0; i < shape[dimension]; i++) {
					_getReflectVector_(strides, shape, dimension, axis,
							vec_pointer, reflected_vec_pointer);
					vec_pointer -= strides[dimension];
				}

				vec_pointer += strides[dimension + 1];

			}

			dimension++;
		}
	}

	void getReflectVector(int*& strides, int*& shape, int dimension, int*& axis,
			double * vec_pointer, double* reflected_vec_pointer) { /*
			 Reflects the array given by pointer vec with strides and shape against axis and write on to reflected_vec.

			 axis must be of the form (0,1,0...,1) etc where 0 means no reflection and 1 means reflection
			 a temporary int initialized to dimension must be supplied.
			 */

		for (int i = 0; i < dimension; i++) {
			if (axis[i] == 1) {
				vec_pointer += (shape[i] - 1) * strides[i];
			}
		}

		if (dimension == 1) {
			if (axis[0] == 0) {
				for (int i = 0; i < shape[0]; i++) {

					*reflected_vec_pointer = *vec_pointer;
					vec_pointer++;
					reflected_vec_pointer++;

				}
			} else {
				for (int i = 0; i < shape[0]; i++) {

					*reflected_vec_pointer = *vec_pointer;
					vec_pointer--;
					reflected_vec_pointer++;

				}
			}

		} else {
			dimension--;
			if (axis[dimension] == 0) {
				for (int i = 0; i < shape[dimension]; i++) {
					_getReflectVector_(strides, shape, dimension, axis,
							vec_pointer, reflected_vec_pointer);
					vec_pointer += strides[dimension];
				}

			} else {

				for (int i = 0; i < shape[dimension]; i++) {
					_getReflectVector_(strides, shape, dimension, axis,
							vec_pointer, reflected_vec_pointer);
					vec_pointer -= strides[dimension];
				}

			}

		}

	}

	void getNextCoordinate(int* array_shape, int* coor, int& dimension) {
		/* Get the next coordinate given by coor for array with array_shape.
		 Write over coor*/

		short i = 0;

		while (i < dimension) {
			if (coor[i] != (array_shape[i] - 1)) {
				coor[i] += 1;
				break;
			} else {
				coor[i] = 0;
				i++;
			}
		}
	}

	void _return_l2_scalar_product_(int& dim, int*& grid_shape,
			int* & gridpt_strides, ELEMENT*& u_p, ELEMENT*& v_p,
			double*& all_gram_values, int* & corner_shifts, int num_corners,
			ELEMENT*&temp_u, ELEMENT*& temp_v, ELEMENT &scalar) {
		/*
		 Procedure:
		 Loop over each grid element
		 Retrieve the corners of u and v at that grid element
		 Calculate the L2 scalar product of u and v at that grid  (by multiplying and summing with the gram values)
		 Add the result to the total scalar product
		 */

		double* agv_p = all_gram_values;
		ELEMENT temp_u_value;

		if (dim == 1) {
			for (int i = 0; i < grid_shape[0]; i++) {
				get_corners(num_corners, u_p, corner_shifts, temp_u);
				get_corners(num_corners, v_p, corner_shifts, temp_v);

				for (int j = 0; j < num_corners; j++) {
					temp_u_value = temp_u[j];

					for (int k = 0; k < num_corners; k++) {
						scalar += temp_u_value * temp_v[k] * *agv_p++;
					}

				}
				u_p++;
				v_p++;
				agv_p = all_gram_values;
			}

			u_p++;
			v_p++;
		} else {

			dim--;

			for (int i = 0; i < grid_shape[dim]; i++)
				_return_l2_scalar_product_(dim, grid_shape, gridpt_strides, u_p,
						v_p, all_gram_values, corner_shifts, num_corners,
						temp_u, temp_v, scalar);

			u_p += gridpt_strides[dim];
			v_p += gridpt_strides[dim];

			dim++;
		}

	}

	void get_grid_info(int* levels, int size, int* grid_shape,
			int* grid_strides, int num_grids, int* gridpt_shape,
			int* gridpt_strides, int num_gridpts) {
		num_grids = 1;
		num_gridpts = 1;

		grid_shape[0] = pow(2, levels[0]);
		gridpt_shape[0] = grid_shape[0] + 1;
		grid_strides[0] = 1;
		gridpt_strides[0] = 1;

		num_grids *= grid_shape[0];
		num_gridpts *= gridpt_shape[0];

		for (int i = 1; i < size; i++) {
			grid_shape[i] = pow(2, levels[i]);
			gridpt_shape[i] = grid_shape[i] + 1;

			grid_strides[i] = grid_strides[i - 1] * grid_shape[i - 1];
			gridpt_strides[i] = gridpt_strides[i - 1] * gridpt_shape[i - 1];

			num_grids *= grid_shape[i];
			num_gridpts *= gridpt_shape[i];

		}
	}

	void _get_corner_shifts_(int* shape, int* strides, int& dimension,
			int& index, int* corner_shifts, int& corner) {

		if (dimension == 1) {

			corner_shifts[corner] = index;
			index++;
			corner++;

			corner_shifts[corner] = index;
			index++;
			corner++;

			index += (shape[0] - 2);
		} else {
			dimension--;

			_get_corner_shifts_(shape, strides, dimension, index, corner_shifts,
					corner);
			_get_corner_shifts_(shape, strides, dimension, index, corner_shifts,
					corner);

			index += (shape[dimension] - 2) * strides[dimension];

			dimension++;
		}
	}

	void get_corner_shifts(int* shape, int* strides, int dimension,
			int* corner_shifts) {
		/*please initialize index, corner=0, reset corner_shifts pointer later too*/

		if (dimension == 1) {
			corner_shifts[0] = 0;
			corner_shifts[1] = 1;

		} else {
			int index = 0;
			int corner = 0;

			dimension--;

			_get_corner_shifts_(shape, strides, dimension, index, corner_shifts,
					corner);
			_get_corner_shifts_(shape, strides, dimension, index, corner_shifts,
					corner);

			dimension++;
		}
	}

	void get_corners(int num_corners, ELEMENT*& u, int*& corner_shifts,
			ELEMENT*& u_corners) {
		/* Having fixed a square grid element in an uniform (possibly isotropic) array of grid values and if the pointer u points to the first corner of that grid element, retrieve the grid values in the corners of the that grid element and put it in u_corners
		 e.g. if 0 1 2 3 4 5 6 7 8 is the array holding

		 0 3 6
		 1 4 7
		 2 5 8

		 and if u points to 4, then the function puts 4 5 7 8 into u_corners.
		 */
		for (int i = 0; i < num_corners; i++) {
			u_corners[i] = u[corner_shifts[i]];
		}

	}

	void copy(double* vec1, double* vec2, int size) {
		// copy values of vec2 to vec1

		for (int i = 0; i < size; i++, vec1++, vec2++)
			*vec1 = *vec2;

	}

	vector<double> getGramValues() {
		vector<double> ret(numPts_, 0.0);
		for (int i = 0; i < numPts_; ++i) {
			ret[i] = gramValues[i];
		}
		return ret;
	}

private:

	int dim_;
	bool * boundaryPoints_;
	int* levels_;
	double* gramValues;
	int numPts_;
};

} /* namespace combigrid */
#endif /* L2SCALARPRODUCT_HPP_ */
