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


struct Parameters{

 /**
  * The perplexity to use in case tsne is the selected algorithm
  */
 int perplexity = 30;

 /**
  * The theta parameter to use in case tsne is the selected algorithm
  */
 float theta=0.5;

 /*
  * The random seed to initialize the selected algorithm
  */
 int seed = 100;

 /*
  * The maximum number of iteration to run on the gradient descent of a selected
  * algorithm
  */
 int maxNunmberIterations = 1000;

 /*
  * The dimentionality to which we want to reduce the data for visualization
  * purposes
  */
 int targetDimension=2;

 /*
  * Number of cores to run the algorithm in parallel
  */
 int numberCores;

};
}

}
