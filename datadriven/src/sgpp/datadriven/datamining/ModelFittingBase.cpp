#include <sgpp/datadriven/datamining/ModelFittingBase.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/application_exception.hpp>


using namespace SGPP::base;

namespace SGPP {
  namespace datadriven {

    ModelFittingBase::ModelFittingBase(SGPP::datadriven::SampleProvider& sampleProvider,
        DataMiningConfiguration config) : sampleProvider(sampleProvider), grid(nullptr), alpha(0) {
    }

    ModelFittingBase::~ModelFittingBase() {
    }

    SGPP::float_t ModelFittingBase::evaluate(DataVector& sample) {
      return 0.0;
    }

    void ModelFittingBase::evaluate(DataMatrix& samples, DataVector& result) {
    }

    OperationMatrix* ModelFittingBase::getRegularizationMatrix(SGPP::datadriven::RegularizationType regType) {
      OperationMatrix* C = NULL;

      if (regType == SGPP::datadriven::RegularizationType::Identity) {
        C = SGPP::op_factory::createOperationIdentity(*grid);
      } else if (regType == SGPP::datadriven::RegularizationType::Laplace) {
        C = SGPP::op_factory::createOperationLaplace(*grid);
      } else {
        throw base::application_exception("ModelFittingBase::getRegularizationMatrix - unknown regularization type");
      }

      return C;
    }

    std::shared_ptr<SGPP::base::Grid> ModelFittingBase::getGrid() {
      return grid;
    }

    std::shared_ptr<SGPP::base::DataVector> ModelFittingBase::getSurpluses() {
      return alpha;
    }

  } /* namespace datadriven */
} /* namespace SGPP */
