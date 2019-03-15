/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */

#ifndef SGPP_PDEFITTER_H
#define SGPP_PDEFITTER_H
#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <list>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
    namespace datadriven {
    class PDFFitter : public ModelFittingDensityEstimationOnOff {
        public:
        std::string configfile;
        explicit PDFFitter(const FitterConfigurationDensityEstimation& conf){
            ModelFittingDensityEstimationOnOff();
            this->config = std::unique_ptr<FitterConfiguration>(
                    std::make_unique<FitterConfigurationDensityEstimation>(conf));
        }
        PDFFitter(){ModelFittingDensityEstimationOnOff();}
        void fit(Dataset &dataset, std::unique_ptr<sgpp::base::Grid> &grid2, std::vector<size_t> ind, bool val);
        };
    }
}

#endif //SGPP_PDEFITTER_H
