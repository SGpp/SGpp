/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Visualizer.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <cstring>

namespace sgpp {

namespace datadriven {

   const VisualizerConfiguration &Visualizer::getVisualizerConfiguration() const{
    return config;
   }


   void Visualizer::createOutputDirectory(unsigned int iteration){

        config.getGeneralConfig().targetFile = config.getGeneralConfig().targetFile+"/Iteration_"
          +std::to_string(iteration);

         std::cout << "Creating output directory " << config.getGeneralConfig().targetFile
           <<std::endl;
         std::string mkdir("mkdir --parents ");

         mkdir.append(config.getGeneralConfig().targetFile);
         system(mkdir.data());

     }
}//namespace datadrivem
}//namespace sgpp
