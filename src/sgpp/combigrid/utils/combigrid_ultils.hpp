/*
 * combigrid_ultils.hpp
 *
 *  Created on: Feb 21, 2011
 *      Author: benk
 */

#ifndef COMBIGRID_ULTILS_HPP_
#define COMBIGRID_ULTILS_HPP_

/** In this hpp file we define the utilities for the combination technique */
#include <math.h>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <complex>
#include <assert.h>

//vector declarations only for internal use , to make the interface more simpler
typedef std::vector<double> DVector;

typedef std::vector<int> IVector;



#define COMBIGRID_OUT(str) { \
	    std::cout << str << std::endl; };

#define COMBIGRID_ERROR_TEST(test,str) { \
	    if (!(test)) std::cout << std::endl << "ERROR: " << str << std::endl; \
	    assert(test); } \

#define COMBIGRID_ERROR_TEST_EQUAL(val1,val2,tolerance,str) { \
	    if (fabs(val1-val2) > tolerance ) {std::cout << std::endl << "ERROR: " << str << \
                   " , v1:" << val1 << " , v2:" << val2 << " , tol:" << tolerance <<std::endl; \
	    assert(false); }} \

#define COMBIGRID_ERROR_EXIT(str) { \
	    std::cout << std::endl << "ERROR: " << str << std::endl; \
	    assert(false); } \

/** write error message out*/
#define COMBIGRID_ERROR_MSG(str) { \
	    std::cout << std::endl << "ERROR: " << str << std::endl; } \

/** Conditionally write out */
#define COMBIGRID_OUT_LEVEL(limit,level,str) \
	   if (limit < level ) COMBIGRID_OUT(str);

/** Conditionally write out */
#define COMBIGRID_OUT_LEVEL1(level,str) \
	   if (level > 1) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL2(level,str) \
	   if (level > 2) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL3(level,str) \
	   if (level > 3) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL4(level,str) \
	   if (level > 4) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL5(level,str) \
	   if (level > 5) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL6(level,str) \
	   if (level > 6) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL7(level,str) \
	   if (level > 7) COMBIGRID_OUT(str);



namespace combigrid{

    /** vector with power two */
	const int powerOfTwo[30] = { 1 , 2 , 4 , 8 , 16 , 32 , 64 , 128 , 256 , 512 , 1024 , 2048 , 4096 ,
			8192 , 16384 , 32768 , 65536 , 131072 , 262144 , 524288 , 1048576 , 2097152 ,
			4194304 , 8388608 , 16777216 , 33554432 , 67108864 , 134217728 , 268435456 };

    /** vector with one over power two */
	const double oneOverPowOfTwo[30]= { 1.0/1.0 , 1.0/2.0 , 1.0/4.0 , 1.0/8.0 , 1.0/16.0 , 1.0/32.0 ,
			1.0/64.0 , 1.0/128.0 , 1.0/256.0 , 1.0/512.0 , 1.0/1024.0 , 1.0/2048.0 , 1.0/4096.0 ,
			1.0/8192.0 , 1.0/16384.0 , 1.0/32768.0 , 1.0/65536.0 , 1.0/131072.0 , 1.0/262144.0 ,
			1.0/524288.0 , 1.0/1048576.0 , 1.0/2097152.0 , 1.0/4194304.0 , 1.0/8388608.0 , 1.0/16777216.0 ,
			1.0/33554432.0 , 1.0/67108864.0 , 1.0/134217728.0 , 1.0/268435456.0 };


	/** the function C_{N}^K , combination of N,K*/
	static int combination(int n, int k)
	{
		if ((k==0)||(n==k)) return 1;
		else
		if ((k==1)||(n==k+1)) return n;
		else
		if ((k==2)||(n==k+2)) return n*(n-1)/2;
		else
		return combination(n-1,k)+combination(n-1,k-1);
	}


}

#endif /* COMBIGRID_ULTILS_HPP_ */
