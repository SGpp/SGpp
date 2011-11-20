/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TEST_dataset_HPP
#define TEST_dataset_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/algorithm/GetAffectedBasisFunctions.hpp"

#include <vector>
#include <utility>
#include <iostream>

namespace sg {
namespace datadriven {

  /**
   * Returns the number of correctly classified instances in data without boundaries
   *
   * @param storage base::GridStorage object that contains the grid points
   * @param basis reference to class that implements to current basis
   * @param alpha the coefficients of the grid points
   * @param data the data the should be tested
   * @param classes the reference classes
   */
  template<class BASIS>
  double test_dataset( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes)
  {
    typedef std::vector<std::pair<size_t, double> > IndexValVector;

    double correct = 0;

#pragma omp parallel shared(correct)
    {
      size_t size = data.getNrows();

      std::vector<double> point;

      base::GetAffectedBasisFunctions<BASIS> ga(storage);

#pragma omp for schedule(static)
      for(size_t i = 0; i < size; i++)
        {

          IndexValVector vec;
          double result = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
            {
              result += iter->second * alpha[iter->first];
            }

          if( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) )
            {
#pragma omp critical
              {
                correct++;
              }
            }
        }
    }

    return correct;
  }

  /**
   * Returns the MSE for given functions values at the evaluation points
   *
   * @param storage base::GridStorage object that contains the grid points
   * @param basis reference to class that implements to current basis
   * @param alpha the coefficients of the grid points
   * @param data the data the should be tested
   * @param refValues the function values at the evaluation points
   */
  template<class BASIS>
  double test_dataset_mse( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues)
  {
    typedef std::vector<std::pair<size_t, double> > IndexValVector;
    base::DataVector result(refValues.getSize());
    double mse = 0;

#pragma omp parallel shared(result)
    {
      size_t size = data.getNrows();
      std::vector<double> point;
      base::GetAffectedBasisFunctions<BASIS> ga(storage);

#pragma omp for schedule(static)
      for(size_t i = 0; i < size; i++)
        {

          IndexValVector vec;
          double res = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
            {
              res += iter->second * alpha[iter->first];
            }

          result[i] = res;
        }
    }

    result.sub(refValues);
    result.sqr();
    mse = result.sum();
    mse /= static_cast<double>(result.getSize());

    return mse;
  }

  /**
   * Returns the number of correctly classified instances in data without boundaries
   *
   * @param storage base::GridStorage object that contains the grid points
   * @param basis reference to class that implements to current basis
   * @param alpha the coefficients of the grid points
   * @param data the data the should be tested
   * @param classes the reference classes
   * @param charaNumbers the number of true positives, true negatives, false positives, false negatives (Vector of length 4)
   * @param ROC_curve DataMatrix into which the ROC curve is stored
   */
  template<class BASIS>
  double test_datasetWithCharacteristicNumber( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers, base::DataMatrix& ROC_curve)
  {
    typedef std::vector<std::pair<size_t, double> > IndexValVector;

    double correct = 0;
    double tp = 0;
    double tn = 0;
    double fp = 0;
    double fn = 0;

#pragma omp parallel shared(correct, tp, tn, fp, fn, ROC_curve)
    {
      size_t size = data.getNrows();

      std::vector<double> point;

      base::GetAffectedBasisFunctions<BASIS> ga(storage);

#pragma omp for schedule(static)
      for(size_t i = 0; i < size; i++)
        {

          IndexValVector vec;
          double result = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
            {
              result += iter->second * alpha[iter->first];
            }

          if( (result >= 0 && classes[i] >= 0) )
            {
#pragma omp critical
              {
                tp++;
                correct++;
              }
            }
          else if( (result < 0 && classes[i] < 0) )
            {
#pragma omp critical
              {
                tn++;
                correct++;
              }
            }
          else if( (result >= 0 && classes[i] < 0) )
            {
#pragma omp critical
              {
                fp++;
              }
            }
          else // ( (result < 0 && classes[i] >= 0) )
            {
#pragma omp critical
              {
                fn++;
              }
            }

          	ROC_curve.set(i, 0, (fp/(fp + tn)));
          	ROC_curve.set(i, 1, (tp/(tp + fn)));
        }
    }

    if (charaNumbers.getSize() < 4)
      {
        charaNumbers.resize(4);
      }

    charaNumbers.set(0, tp);
    charaNumbers.set(1, tn);
    charaNumbers.set(2, fp);
    charaNumbers.set(3, fn);

    return correct;
  }

}
}

#endif /* TEST_dataset_HPP */
