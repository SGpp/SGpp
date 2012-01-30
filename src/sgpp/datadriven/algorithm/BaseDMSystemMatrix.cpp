/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/algorithm/BaseDMSystemMatrix.hpp"

namespace sg
{
namespace datadriven
{

BaseDMSystemMatrix::BaseDMSystemMatrix(sg::base::DataMatrix& trainData, double lambda) : data_(&trainData), lambda_(lambda)
{
}

}

}
