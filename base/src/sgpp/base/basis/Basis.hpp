#ifndef BASIS_HPP
#define BASIS_HPP

#include <sgpp/globaldef.hpp>


namespace SGPP{
namespace base {

	template<class LT, class IT>
	class Basis {
      public:
		virtual double eval(LT level, IT index, double p) = 0;
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
