// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/domain/StretchingFactory.hpp>

namespace combigrid {

AbstractStretchingMaker* createStretchingMaker(combigrid::Stretching str) {
  AbstractStretchingMaker* maker = NULL;

  switch (str) {
    case ATAN:
      maker = new AtanSpecialStretching();
      break;

    case TAN:
      maker = new TanStretching();
      break;

    case CHEBYSHEV:
      maker = new CombiChebyshevStretching();
      break;

    case EQUIDISTANT:
      maker = new CombiEquidistantStretching();
      break;

    case LEGENDRE:
      maker = new CombiLegendreStretching();
      break;

    case BASU:
      maker = new CombiBasuStretching();
      break;

    case UNKNOWN:
    default:
      COMBIGRID_OUT_WRN("Input stretching type unknown! Returning a NULL stretching maker!",
                        __FILE__, __LINE__)
      break;
  }

  return maker;
}
}  // namespace combigrid
