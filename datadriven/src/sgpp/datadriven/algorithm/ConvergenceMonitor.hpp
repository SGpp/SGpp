// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONVERGENCEMONITOR_HPP
#define CONVERGENCEMONITOR_HPP

#include <deque>
#include <sgpp/globaldef.hpp>

/**
 * Description
 * 
 * 
 */
class ConvergenceMonitor {
  public:
    /**
    * Constructor
    */
    ConvergenceMonitor(double pDeclineThreshold,
                       size_t pBufferSize, size_t pMinRefInterval);

    /**
    * Destructor
    */
    virtual ~ConvergenceMonitor();
  
    /**
    * Description
    */
    virtual void pushToBuffer(double currentValidError,
                              double currentTrainError);
    
    /**
    * Description
    */
    virtual bool checkConvergence();

    size_t nextRefCnt;
    size_t minRefInterval;
    std::deque<double> validErrorDeclineBuffer;
    std::deque<double> trainErrorDeclineBuffer;

  protected:

  private:
    /*sgpp::base::DataMatrix& trainData;
    sgpp::base::DataMatrix& validData;
    sgpp::base::DataVector& trainLabels;
    sgpp::base::DataVector& validLabels;*/
    double validErrorSum1;
    double validErrorSum2;
    double trainErrorSum1;
    double trainErrorSum2;
    //double validAvgError1;
    //double validAvgError2;
    //double trainAvgError1;
    //double trainAvgError2;
    //double currentValidError;
    //double currentTrainError;
    double validRatio;
    double trainRatio;
    double declineThreshold;
    size_t bufferSize;
    //std::deque<double> validErrorDeclineBuffer;
    //std::deque<double> trainErrorDeclineBuffer;   
	
};

#endif /* CONVERGENCEMONITOR_HPP */
