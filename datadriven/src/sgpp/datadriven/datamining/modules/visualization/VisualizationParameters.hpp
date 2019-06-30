/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossvalidationConfiguration.hpp
 *
 *  Created on: Jun 8, 2019
 *      Author: Vincent Bautista
 */

#pragma once

namespace sgpp {
namespace datadriven {


struct VisualizationParameters{

 /**
  * The perplexity to use in case tsne is the selected algorithm
  */
 double perplexity = 30;

 /**
  * The theta parameter to use in case tsne is the selected algorithm
  */
 double theta=0.5;

 /*
  * The random seed to initialize the selected algorithm
  */
 size_t seed = 100;

 /*
  * The maximum number of iteration to run on the gradient descent of a selected
  * algorithm
  */
 size_t maxNunmberIterations = 1000;

 /*
  * The dimentionality to which we want to reduce the data for visualization
  * purposes
  */
 size_t targetDimension=2;

 /*
  * Number of cores to run the algorithm in parallel
  */
 size_t numberCores=1;

};
}

}
