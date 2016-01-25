#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

using namespace SGPP::base;
using namespace SGPP::datadriven;

namespace SGPP {
  namespace datadriven {

    void OperationDensityMarginalizeKDE::doMarginalize(size_t mdim,
        GaussianKDE& marginalizedKDE) {
      // dimensionality of new set
      size_t ndim = kde->getDim();

      if (ndim <= 1) {
        throw data_exception(
          "OperationDensityMarginalizeKDE::doMarginalize : Dimensionality of kde is <= 1");
      }

      if (mdim >= ndim) {
        throw data_exception(
          "OperationDensityMarginalizeKDE::doMarginalize : Dimensionality is out of range");
      }

      std::vector<DataVector*> newSamplesVec(ndim - 1);

      // copy all values but the ones in dim
      size_t currentDim = 0;

      for (size_t idim = 0; idim < ndim; idim++) {
        if (idim != mdim) {
          newSamplesVec[currentDim] = kde->getSamples(idim);
          currentDim++;
        }
      }

      // initialize kde with new samples
      marginalizedKDE.initialize(newSamplesVec);
    }

    void OperationDensityMarginalizeKDE::doMarginalize(std::vector<size_t>& mdims,
        GaussianKDE& marginalizedKDE) {
      // dimensionality of original and marginalized set
      size_t ndims = kde->getDim();
      size_t newndims = mdims.size();
      // container for new sample set
      std::vector<DataVector*> newSamplesVec(newndims);
      // copy all enries but the ones in mdim
      size_t currentDim = 0;

      for (size_t idim = 0; idim < ndims; idim++) {
        // search if the current dimension in mdims
        size_t jdim = 0;

        while (jdim < newndims) {
          if (jdim == mdims[jdim]) {
            // do marginalization -> skip samples in current dimension
            break;
          }

          jdim++;
        }

        // dimension found which should not be marginalized
        if (jdim == newndims) {
          newSamplesVec[currentDim] = kde->getSamples(idim);
          currentDim++;
        }
      }

      marginalizedKDE.initialize(newSamplesVec);
    }

    void OperationDensityMarginalizeKDE::margToDimX(size_t mdim,
        GaussianKDE& marginalizedKDE) {
      if (mdim >= kde->getDim()) {
        throw new data_exception(
          "OperationDensityMarginalizeKDE::margToDimX; mdim out of range");
      }

      // dimensionality of new set
      std::vector<DataVector*> newSamplesVec(1);
      newSamplesVec[0] = kde->getSamples(mdim);

      // initialize marginalized kde
      marginalizedKDE.initialize(newSamplesVec);
    }

    void OperationDensityMarginalizeKDE::margToDimXs(std::vector<size_t>& mdims,
        GaussianKDE& marginalizedKDE) {
      // dimensionality of new set
      size_t ndimsNew = mdims.size();
      std::vector<DataVector*> newSamplesVec(ndimsNew);

      // get the subset of the data for marginalized density
      for (size_t idim = 0; idim < ndimsNew; idim++) {
        newSamplesVec[idim] = kde->getSamples(mdims[idim]);
      }

      // initialize kde with new samples
      marginalizedKDE.initialize(newSamplesVec);
    }

  }
}
