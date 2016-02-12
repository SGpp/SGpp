// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STRETCHINGFACTORY_HPP_
#define STRETCHINGFACTORY_HPP_

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>
#include <sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp>
#include <sgpp/combigrid/domain/CombiBasuStretching.hpp>
#include <sgpp/combigrid/domain/CombiChebyshevStretching.hpp>
#include <sgpp/combigrid/domain/CombiEquidistantStretching.hpp>
#include <sgpp/combigrid/domain/CombiLegendreStretching.hpp>
#include <sgpp/combigrid/domain/CombiTanStretching.hpp>

namespace combigrid {

AbstractStretchingMaker* createStretchingMaker(combigrid::Stretching str);
//
//    AbstractStretchingMaker* maker = NULL;
//
//    switch (str){
//    case ATAN:
//      maker = new AtanSpecialStretching();
//      break;
//    case TAN:
//      maker = new TanStretching();
//      break;
//    case CHEBYSHEV:
//      maker = new CombiChebyshevStretching();
//      break;
//    case EQUIDISTANT:
//      maker = new CombiEquidistantStretching();
//      break;
//    case LEGENDRE:
//      maker = new CombiLegendreStretching();
//      break;
//    case BASU:
//      maker = new CombiBasuStretching();
//      break;
//    default:
//      break;
//    }
//
//    return maker;
//  }
}  // namespace combigrid

#endif /* STRETCHINGFACTORY_HPP_ */
