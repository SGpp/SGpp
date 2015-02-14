// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef BASIS_HPP
#define BASIS_HPP

#include <sgpp/globaldef.hpp>


namespace SGPP{
namespace base {

	template<class LT, class IT>
	class Basis {
      public:
		virtual float_t eval(LT level, IT index, float_t p) = 0;
		virtual ~Basis(){};
		//Basis();
      /*private:
		Basis(Basis const&);
		Basis& operator=(Basis const&);
		*/
	};

	typedef Basis<unsigned int, unsigned int> SBasis;
}
}

#endif // BASIS_HPP
