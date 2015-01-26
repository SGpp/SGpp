/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#include <cstdlib>
#include <ctime>
#include "Random.hpp"

#ifdef USE_STD_RANDOM
#include <random>
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace mcm {
      
    bool Random::is_seeded = false;
#ifdef USE_STD_RANDOM
    std::mt19937 Random::gen = std::mt19937();
    std::uniform_int_distribution<int> Random::dist = std::uniform_int_distribution<int>(0, RAND_MAX);
#endif


    void Random::seed(int seed_value) {
#ifdef USE_STD_RANDOM       
        gen.seed(seed_value);
#else
        srand(seed_value);
#endif
        is_seeded = true;
    }
    
    int Random::random() {
        if (!is_seeded) {
            seed((int)std::time(NULL));
        }
#ifdef USE_STD_RANDOM 
        return dist(gen);
#else
        return rand();
#endif   
    }
    
    double Random::random_double() {
        return (double)(random())/RAND_MAX;
    }

  }
}

