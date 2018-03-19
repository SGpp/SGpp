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
  void createRandomConfigs(size_t nBits, std::vector<int>& configIDs, int seed, size_t start);
  void calculateConstrainedSpace(const DataVector& transformedScores, double lambda, int shrink);
  void transformScores(const DataVector& source, DataVector& target);
  bool fixConfigBits();
  void resetBits();
  void setParameters(int configID, int matrixrow);
  void addConstraint(int idx, int bias);
  bool checkConstraints();
  int moveToNewSpace(int configID, std::vector<ConfigurationBit*> oldFreeBits);



  // double expkernel(base::DataVector x1, base::DataVector x2);
protected:
	base::DataMatrix paritymatrix;
  FitterFactory* fitterFactory;
  std::vector<int> configIDs;
  DataVector savedScores;
  std::vector<std::vector<ConfigurationBit*> > parityrow;
  std::vector<ConfigurationBit*> freeBits;
  std::vector<ConfigurationBit*> configBits;
  std::vector<std::unique_ptr<ConfigurationRestriction>> constraints;


};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_ */
