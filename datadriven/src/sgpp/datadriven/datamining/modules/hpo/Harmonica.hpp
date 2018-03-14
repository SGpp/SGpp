/*
 * Harmonica.hpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include "FitterFactory.hpp"


namespace sgpp {
namespace datadriven {

class Harmonica {
public:
	Harmonica(FitterFactory* fitterFactory);


  void prepareConfigs(std::vector<ModelFittingBase*>& fitters);
  void createRandomConfigs(int nBits, std::vector<int>& configIDs, int seed);
  void calculateConstrainedSpace(const DataVector& transformedScores, int lambda, int shrink);
  void transformScores(const DataVector& source, DataVector& target);



    // double expkernel(base::DataVector x1, base::DataVector x2);
protected:
	base::DataMatrix paritymatrix;
  int nBits;
  FitterFactory* fitterFactory;


};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_ */
