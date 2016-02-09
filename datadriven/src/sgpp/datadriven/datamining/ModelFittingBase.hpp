// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#pragma once

#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/datadriven/datamining/SampleProvider.hpp>


namespace SGPP {
  namespace datadriven {

    class ModelFittingBase {
      public:
        ModelFittingBase(SGPP::datadriven::SampleProvider& sampleProvider);

        virtual ~ModelFittingBase();

        /**
         *
         */
        virtual void fit() = 0;

        /**
         *
         * @param sample
         * @return
         */
        virtual SGPP::float_t evaluate(SGPP::base::DataVector& sample);

        /**
         *
         * @param samples
         * @return
         */
        virtual void evaluate(SGPP::base::DataMatrix& samples, SGPP::base::DataVector& result);

        virtual std::shared_ptr<SGPP::base::Grid> getGrid();
        virtual std::shared_ptr<SGPP::base::DataVector> getSurpluses();

      protected:
        SGPP::base::OperationMatrix* getRegularizationMatrix(SGPP::datadriven::RegularizationType regType);

        SGPP::datadriven::SampleProvider& sampleProvider;
        std::shared_ptr<SGPP::base::Grid> grid;
        std::shared_ptr<SGPP::base::DataVector> alpha;
    };

  } /* namespace datadriven */
} /* namespace SGPP */

