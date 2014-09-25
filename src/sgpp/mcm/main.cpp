#include "sgpp_mcm.hpp"
#include <iostream>

using namespace sg::base;
using namespace sg::mcm;

    int main(void){
      size_t dim = 1;
      
      NaiveSampleGenerator* nsg = new NaiveSampleGenerator(dim);
      DataVector dv1(dim);
      nsg->getSample(dv1);
      
      std::cout << dv1.toString() << std::endl;
      
      /* Stratified Sampling Test */
      
      const size_t dim_ss = 2;
      /*
      long long int unterteilungen[dim_ss] = { 3 , 3 };
      StratifiedSampleGenerator* ssg = new StratifiedSampleGenerator( dim_ss, unterteilungen );
      
      DataVector dv2( dim_ss );
      
      std::cout << "Stratified Samples" << std::endl;
      std::cout << "Anzahl Samples: " << ssg->getNumberOfSamples() << std::endl;
      
      while( ssg->hasNext() ) {
	ssg->getNext( dv2 );
	std::cout << dv2.toString() << std::endl;
      }
      
      /* Latin Hypercube Sampling Test */
      /*
      std::cout << "Latin Hypercube Samples" << std::endl;
      
      const size_t dim_lhs = 10;
      LatinHypercubeSampleGenerator* lhsg = new LatinHypercubeSampleGenerator( dim_lhs, 10 );
      
      DataVector dv3( dim_lhs );
      
      std::cout << "Anzahl Samples: " << lhsg->getNumberOfSamples() << std::endl;
      
      while( lhsg->hasNext() ) {
	lhsg->getNext( dv3 );
	std::cout << dv3.toString() << std::endl;
      }*/
      
      std::cout << "Quasi Sobol Sequence Scrambled Samples" << std::endl;
      
      /*size_t seed = 0;
      SobolSequenceGenerator* qsg = new SobolSampleGenerator(dim,seed);
      DataVector dv_q(dim);
      qsg->getSample( dv_q );
      /*for(size_t i = 0; i<1; i++ ) {
	qsg->getSample( dv_q );
	std::cout << dv_q.toString() << std::endl;
      }*/
      return(0);
      
    }
