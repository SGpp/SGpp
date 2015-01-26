// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef LEARNERVECTORIZEDIDENTITYSP_HPP
#define LEARNERVECTORIZEDIDENTITYSP_HPP

#include <sgpp/datadriven/application/LearnerBaseSP.hpp>

#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace parallel {

    /**
     * This class implements standard sparse grid regression
     * with an Identity matrix as regularization operator.
     *
     * Furthermore this Learner provides support for several
     * vectorization approaches covering GPUs, CPUs and coprocessors.
     *
     * This version supports single precision floating point numbers.
     */
    class LearnerVectorizedIdentitySP : public SGPP::datadriven::LearnerBaseSP {
      protected:
        /// vectorization selector
        VectorizationType vecType_;

        MPIType mpiType_;

        virtual SGPP::datadriven::DMSystemMatrixBaseSP* createDMSystem(SGPP::base::DataMatrixSP& trainDataset, float lambda);

        virtual void postProcessing(const SGPP::base::DataMatrixSP& trainDataset, const SGPP::solver::SLESolverType& solver,
                                    const size_t numNeededIterations);

      public:
        /**
         * Constructor
         *
         * @param vecType selection of vectorization to employ
         * @param isRegression set to true if a regression task should be executed
         * @param isVerbose set to true in order to allow console output
         */
        LearnerVectorizedIdentitySP(const VectorizationType vecType, const bool isRegression, const bool isVerbose = true);
        LearnerVectorizedIdentitySP(const VectorizationType vecType, const MPIType mpiType, const bool isRegression, const bool isVerbose = true);

        /**
         * Constructor
         *
         * @param tGridFilename path to file that contains a serialized grid
         * @param tAlphaFilename path to file that contains the grid's coefficients
         * @param vecType selection of vectorization to employ
         * @param isRegression set to true if a regression task should be executed
         * @param isVerbose set to true in order to allow console output
         */
        LearnerVectorizedIdentitySP(const std::string tGridFilename, const std::string tAlphaFilename, const VectorizationType vecType,
                                    const bool isRegression, const bool isVerbose = true);

        /**
         * Destructor
         */
        virtual ~LearnerVectorizedIdentitySP();

        virtual SGPP::base::DataVectorSP predict(SGPP::base::DataMatrixSP& testDataset);
    };

  }

}

#endif /* LEARNERVECTORIZEDIDENTITYSP_HPP */