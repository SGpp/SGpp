// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "Random.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "ScrambledSobolSampleGenerator.hpp"

using namespace SGPP::base;
using namespace std;

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace mcm {

ScrambledSobolSampleGenerator::ScrambledSobolSampleGenerator(
		size_t dimensions) :
		SampleGenerator(dimensions) {
	this->seed = 1;
	this->numberOfCurrentSample = 1;
}

ScrambledSobolSampleGenerator::ScrambledSobolSampleGenerator(
		size_t dimensions, size_t seed) :
		SampleGenerator(dimensions) {
	this->seed = seed + 1;
	this->numberOfCurrentSample = 1;
}

size_t ScrambledSobolSampleGenerator::i4_bit_hi1(size_t n) {
	double i = floor(n);
	size_t bit = 0;
	while (true) {
		if (i <= 0)
			break;
		bit += 1;
		i = floor(i / 2.);
	}
	return bit;
}

size_t ScrambledSobolSampleGenerator::i4_bit_lo0(size_t n) {
	size_t bit = 0;
	double i = floor(n);
	double i2;
	while (true) {
		bit = bit + 1;
		i2 = floor(i / 2.);
		if (i == 2 * i2)
			break;
		i = i2;
	}
	return bit;
}

size_t dim_num_save_scrambled;
bool initialized_scrambled = false;
size_t seed_save_scrambled;
size_t dim_max_scrambled;
size_t log_max_scrambled;
size_t atmost_scrambled;
std::vector<size_t> lastq_scrambled;
size_t maxcol_scrambled;
size_t poly_scrambled[] = { 1, 3, 7, 11, 13, 19, 25, 37, 59, 47, 61, 55, 41, 67,
		97, 91, 109, 103, 115, 131, 193, 137, 145, 143, 241, 157, 185, 167, 229,
		171, 213, 191, 253, 203, 211, 239, 247, 285, 369, 299 };
double recipd_scrambled;
SGPP::base::DataMatrix* v_scrambled;

void ScrambledSobolSampleGenerator::getSample(DataVector& dv) {

	size_t dim_num = dv.getSize();

	if (!initialized_scrambled or dim_num != dim_num_save_scrambled) {
		initialized_scrambled = true;
		dim_max_scrambled = 40;
		dim_num_save_scrambled = -1;
		log_max_scrambled = 30;
		seed_save_scrambled = -1;
		//
		//    Initialize (part of) V.
		//
		v_scrambled = new SGPP::base::DataMatrix((size_t) dim_max_scrambled,
				log_max_scrambled);

		for (size_t i = 0; i < dim_max_scrambled; i++) {
			for (size_t j = 0; j < log_max_scrambled; j++) {
				v_scrambled->SGPP::base::DataMatrix::set(i, j, 0);
			}
		}

		//v_scrambled[0:40,0]
		double d0[] = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
				1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
				1., 1., 1., 1., 1., 1., 1., 1., 1., 1. };
		DataVector temp_dv0(d0, 40);
		v_scrambled->setColumn(0, temp_dv0);

		//v[2:40,1]
		double d1[] = { 0, 0, 1, 3, 1, 3, 1, 3, 3, 1, 3, 1, 3, 1, 3, 1, 1, 3, 1,
				3, 1, 3, 1, 3, 3, 1, 3, 1, 3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 3 };
		DataVector temp_dv_scrambled(d1, 40);
		v_scrambled->setColumn(1, temp_dv_scrambled);

		//v[3:40,2]
		double d2[] = { 0, 0, 0, 7, 5, 1, 3, 3, 7, 5, 5, 7, 7, 1, 3, 3, 7, 5, 1,
				1, 5, 3, 3, 1, 7, 5, 1, 3, 3, 7, 5, 1, 1, 5, 7, 7, 5, 1, 3, 3 };
		DataVector temp_dv2(d2, 40);
		v_scrambled->setColumn(2, temp_dv2);

		//v[5:40,3]
		double d3[] = { 0, 0, 0, 0, 0, 1, 7, 9, 13, 11, 1, 3, 7, 9, 5, 13, 13,
				11, 3, 15, 5, 3, 15, 7, 9, 13, 9, 1, 11, 7, 5, 15, 1, 15, 11, 5,
				3, 1, 7, 9 };
		DataVector temp_dv3(d3, 40);
		v_scrambled->setColumn(3, temp_dv3);

		//v[7:40,4]
		double d4[] = { 0, 0, 0, 0, 0, 0, 0, 9, 3, 27, 15, 29, 21, 23, 19, 11,
				25, 7, 13, 17, 1, 25, 29, 3, 31, 11, 5, 23, 27, 19, 21, 5, 1,
				17, 13, 7, 15, 9, 31, 9 };
		DataVector temp_dv4(d4, 40);
		v_scrambled->setColumn(4, temp_dv4);

		//v[13:40,5]
		double d5[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 37, 33, 7, 5, 11,
				39, 63, 27, 17, 15, 23, 29, 3, 21, 13, 31, 25, 9, 49, 33, 19,
				29, 11, 19, 27, 15, 25 };
		DataVector temp_dv5(d5, 40);
		v_scrambled->setColumn(5, temp_dv5);

		//v[19:40,6]
		double d6[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				13, 33, 115, 41, 79, 17, 29, 119, 75, 73, 105, 7, 59, 65, 21, 3,
				113, 61, 89, 45, 107 };
		DataVector temp_dv6(d6, 40);
		v_scrambled->setColumn(6, temp_dv6);

		//v[37:40,7]
		double d7[] =
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
						0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 23,
						39 };
		DataVector temp_dv7(d7, 40);
		v_scrambled->setColumn(7, temp_dv7);

		atmost_scrambled = (size_t) (pow(2., (double) log_max_scrambled) - 1);
		//
		//    Find the number of bits in ATMOST.
		//
		maxcol_scrambled = i4_bit_hi1(atmost_scrambled);

		//
		//    Initialize row 1 of V.
		//
		for (size_t i = 0; i <= maxcol_scrambled; i++)
			v_scrambled->set(0, i, 1);
	}
	//
	//    Things to do only if the dimension changed.
	//

	if (dim_num != dim_num_save_scrambled) {
		//
		//    Check parameters.
		//
		if (dim_num < 1 or dim_max_scrambled < dim_num) {
			std::cout << "I4_SOBOL - Fatal error!" << std::endl;
			std::cout << "The spatial dimension DIM_NUM should satisfy:"
					<< std::endl;
			std::cout << "1 <= DIM_NUM <= dim_max1" << std::endl;
			std::cout << "But this input value is DIM_NUM ..." << std::endl;
			return;
		}

		dim_num_save_scrambled = dim_num;
		//
		//    Initialize the remaining rows of V.
		//
		for (size_t i = 2; i < dim_num + 1; i++) {
			//
			//    The bits of the integer POLY(I) gives the form of poly1nomial I.
			//
			//    Find the degree of poly1nomial I from binary encoding.
			//
			size_t j = poly_scrambled[i - 1];
			size_t m = 0;
			while (true) {
				j = (size_t) floor(j / (size_t) 2);
				if (j <= 0)
					break;
				m = m + 1;
			}
			//
			//    Expand this bit pattern to separate components of the logical array INCLUD.
			//
			j = poly_scrambled[i - 1];
			std::vector<size_t> includ(m);
			size_t j2;
			for (size_t k = m; k > 0; k--) {		//k um +1 geshifted
				j2 = (size_t) floor(j / (size_t) 2);
				includ[k - 1] = (j != 2 * j2);
				j = j2;
			}

			//
			//    Calculate the remaining elements of row I as explained
			//    in Bratley and Fox, section 2.
			//

			size_t newv;
			size_t l;
			for (size_t j = m + 1; j < maxcol_scrambled + 1; j++) {
				newv = (size_t) v_scrambled->get(i - 1, j - m - 1);
				l = 1;
				for (size_t k = 1; k < m + 1; k++) {
					l = 2 * l;
					if (includ[k - 1]) {
						newv = newv
								^ (l
										* (size_t) v_scrambled->get(i - 1,
												j - k - 1));
					}
				}
				v_scrambled->set(i - 1, j - 1, (double) newv);
			}
		}

		//
		//    Multiply columns of V by appropriate power of 2.
		//

		double l = 1;
		for (size_t j = maxcol_scrambled - 1; j > 0; j--) {
			l = 2 * l;
			for (size_t i = 0; i < dim_num; i++) {
				v_scrambled->set(i, j - 1, v_scrambled->get(i, j - 1) * l);
			}
		}

		//
		//    RECIPD is 1/(common denominator of the elements in V).
		//
		recipd_scrambled = 1.0 / (2.0 * l);
		lastq_scrambled.reserve(dim_num);
		for (size_t i = 0; i < dim_num; i++)
			lastq_scrambled[i] = 0;
	}

	seed = (size_t) floor(seed);
	size_t l = 0;
	if (seed < 0)
		seed = 0;

	if (seed == 0) {
		l = 1;
		lastq_scrambled.reserve(dim_num);
		for (size_t i = 0; i < dim_num; i++)
			lastq_scrambled[i] = 0;
	} else if (seed == seed_save_scrambled + 1) {
		//
		//    Find the position of the right-hand zero in SEED.
		//
		l = i4_bit_lo0(seed);

	} else if (seed <= seed_save_scrambled) {
		seed_save_scrambled = 0;
		l = 1;
		lastq_scrambled.reserve(dim_num);
		for (size_t i = 0; i < dim_num; i++)
			lastq_scrambled[i] = 0;

		for (size_t seed_temp = seed_save_scrambled; seed_temp < seed;
				seed_temp++) {
			l = i4_bit_lo0(seed_temp);
			for (size_t i = 1; i < dim_num + 1; i++) {
				lastq_scrambled[i - 1] = lastq_scrambled[i - 1]
						^ (size_t) v_scrambled->get(i - 1, l - 1);
			}
		}
		l = i4_bit_lo0(seed);
	} else if (seed_save_scrambled + 1 < seed) {
		for (size_t seed_temp = seed_save_scrambled + 1; seed_temp < seed;
				seed_temp++) {
			l = i4_bit_lo0(seed_temp);
			for (size_t i = 1; i < dim_num + 1; i++) {
				lastq_scrambled[i - 1] = lastq_scrambled[i - 1]
						^ (size_t) v_scrambled->get(i - 1, l - 1);
			}
		}
		l = i4_bit_lo0(seed);
	}

	//
	//    Check that the user is not calling too many times!
	//
	if (maxcol_scrambled < l) {
		std::cout << "I4_SOBOL - Fatal error!" << std::endl;
		std::cout << "Too many calls!" << std::endl;
		std::cout << "MAXCOL = %d\n'%maxcol_scrambled" << std::endl;
		std::cout << "L =            %d\n'%l" << std::endl;
	}
	//
	//    Calculate the new components of QUASI.
	//

	for (size_t i = 1; i < dim_num + 1; i++) {
		dv[i - 1] = (double) lastq_scrambled[i - 1] * recipd_scrambled;
		lastq_scrambled[i - 1] = lastq_scrambled[i - 1]
				^ (size_t) v_scrambled->get(i - 1, l - 1);
	}

	seed_save_scrambled = seed;
	seed = seed + 1;
	scramble(dv);
}

void ScrambledSobolSampleGenerator::permutation(int base, int* permutatuion) {
	// Permutation bestimmen
	for (size_t temp = 0; (int) temp < base; temp++)
		permutatuion[temp] = -1;

	int steps = 0;
	int position = 0;

	for (int i = base - 1; i >= 0; i--) {
		steps = Random::random() % base;  //0-9, gleichverteilt
		while (steps >= 0) {
			steps--;
			if (permutatuion[(position + 1) % base] > -1) {
				while (permutatuion[position] > -1) {
					position = (position + 1) % base;
				}
			} else {
				position = (position + 1) % base;
			}
		}
		permutatuion[position] = i;
	}
}

void ScrambledSobolSampleGenerator::scramble(DataVector& dv) {

	//Permutieren
	int base = 8;
	int* permutatuion = new int[base];
	permutation(base, permutatuion);

	//Extrahiere die einzelnen Zeichen
	for (size_t i = 0; i < dv.getSize(); i++) {
		double zahl = dv[i];
		double exp = pow((double) base, -1.0);		//Basis 8

		double a = zahl / exp;
		double b = floor(zahl / exp);
		double exp_inc = exp;

		double ausgabezahl = 0.0;
		size_t j = 0;
		while (((a - b) > 0.000001) && j < 13) {//End after calculating 13 digits
			j = j + 1;

			ausgabezahl = ausgabezahl
					+ ((double) permutatuion[(size_t) floor(a)] * exp_inc);
			exp_inc = exp_inc * exp;

			a = a - b;
			a = a / exp;
			b = floor(a);
		}
		dv[i] = ausgabezahl;

		//dv[i] = 0.5;

	}
}

}
}