/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TEST_dataset_HPP
#define TEST_dataset_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>

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
    double test_dataset( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      double correct = 0;

      #pragma omp parallel shared(correct)
      {
        size_t size = data.getNrows();

        std::vector<double> point;

        base::GetAffectedBasisFunctions<BASIS> ga(storage);

        #pragma omp for schedule(static)

        for (size_t i = 0; i < size; i++) {

          IndexValVector vec;
          double result = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
            result += iter->second * alpha[iter->first];
          }

          if ( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) ) {
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
    double test_dataset_mse( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;
      base::DataVector result(refValues.getSize());
      double mse = 0;

      #pragma omp parallel shared(result)
      {
        size_t size = data.getNrows();
        std::vector<double> point;
        base::GetAffectedBasisFunctions<BASIS> ga(storage);

        #pragma omp for schedule(static)

        for (size_t i = 0; i < size; i++) {

          IndexValVector vec;
          double res = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
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
     * @param threshold threshold which decides if an instance belongs a given class
     */
    template<class BASIS>
    double test_datasetWithCharacteristicNumber( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers, double threshold) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      double correct = 0;
      double tp = 0;
      double tn = 0;
      double fp = 0;
      double fn = 0;

      #pragma omp parallel shared(correct, tp, tn, fp, fn)
      {
        size_t size = data.getNrows();

        std::vector<double> point;

        base::GetAffectedBasisFunctions<BASIS> ga(storage);

        #pragma omp for schedule(static)

        for (size_t i = 0; i < size; i++) {

          IndexValVector vec;
          double result = 0;

          data.getRow(i, point);

          ga(basis, point, vec);

          for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
            result += iter->second * alpha[iter->first];
          }

          if ( (result >= threshold && classes[i] >= 0) ) {
            #pragma omp critical
            {
              tp++;
              correct++;
            }
          } else if ( (result < threshold && classes[i] < 0) ) {
            #pragma omp critical
            {
              tn++;
              correct++;
            }
          } else if ( (result >= threshold && classes[i] < 0) ) {
            #pragma omp critical
            {
              fp++;
            }
          } else { // ( (result < threshold && classes[i] >= 0) )
            #pragma omp critical
            {
              fn++;
            }
          }
        }
      }

      if (charaNumbers.getSize() < 4) {
        charaNumbers.resize(4);
      }

      charaNumbers.set(0, tp);
      charaNumbers.set(1, tn);
      charaNumbers.set(2, fp);
      charaNumbers.set(3, fn);

      return correct;
    }

    /**
     * Returns the number of correctly classified instances in data without boundaries
     *
     * @param storage base::GridStorage object that contains the grid points
     * @param basis reference to class that implements to current basis
     * @param alpha the coefficients of the grid points
     * @param data the data the should be tested
     * @param classes the reference classes
     * @param thresholds the thresholds (between -1.0 and 1.0) for calculating the ROC curve
     * @param ROC_curve DataMatrix into which the ROC curve is stored
     */
    template<class BASIS>
    void test_calculateROCcurve( base::GridStorage* storage, BASIS& basis, base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) {
      size_t num_points = thresholds.getSize();

      if (ROC_curve.getNrows() != num_points) {
        ROC_curve.resize(num_points);
      }

      base::DataVector charNum(4);

      for (size_t i = 0; i < num_points; i++) {
        test_datasetWithCharacteristicNumber( storage, basis, alpha, data, classes, charNum, thresholds.get(i));

        double tp = charNum.get(0);
        double tn = charNum.get(1);
        double fp = charNum.get(2);
        double fn = charNum.get(3);

        // 1-spec.
        ROC_curve.set(i, 0, (fp / (fp + tn)));
        // sensi.
        ROC_curve.set(i, 1, (tp / (tp + fn)));
      }
    }
  }
}

#endif /* TEST_dataset_HPP */
