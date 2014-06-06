/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/common/ArBBKernels5D.hpp"

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <string>

#include <arbb.hpp>

namespace sg {

  namespace parallel {

    typedef arbb::array<arbb::f64, 5> vecElem64;
    typedef arbb::uncaptured< arbb::array<double, 5> >::type uc_vecElem64;
    typedef arbb::array<arbb::f32, 5> vecElem32;
    typedef arbb::uncaptured< arbb::array<float, 5> >::type uc_vecElem32;

    arbb::dense< vecElem32, 1> ArBB_DataSP_5D;
    arbb::dense< vecElem32, 1> ArBB_LevelSP_5D;
    arbb::dense< vecElem32, 1> ArBB_IndexSP_5D;

    arbb::dense< vecElem64, 1> ArBB_Data_5D;
    arbb::dense< vecElem64, 1> ArBB_Level_5D;
    arbb::dense< vecElem64, 1> ArBB_Index_5D;

    template <typename fp_Type>
    void arbb_multTrans(const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPointTrans(const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Data, const arbb::array<fp_Type, 5>& Level, const arbb::array<fp_Type, 5>& Index, const arbb::dense<fp_Type>& source, fp_Type& result) {
          result = 0.0;

          _for (arbb::usize i = 0, i < Data.num_cols(), i++) {
            result += (source[i] * arbb::mul_reduce(arbb::max((1.0 - arbb::abs(((Level * Data[i]) - Index))), 0.0)));
          }
          _end_for;
        }
      };

      arbb::map(local::evalGridPointTrans)(Data, Level, Index, source, result);
    }

    template <typename fp_Type>
    void arbb_mult(const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPoint(const arbb::array<fp_Type, 5>& DataPoint, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 5>, 1 >& Index, const arbb::dense<fp_Type>& alpha, fp_Type& result) {
          result = 0.0;

          _for (arbb::usize j = 0, j < Level.num_cols(), j++) {
            result += (alpha[j] * arbb::mul_reduce(arbb::max((1.0 - arbb::abs(((Level[j] * DataPoint) - Index[j]))), 0.0)));
          }
          _end_for;
        }
      };

      arbb::map(local::evalGridPoint)(Data, Level, Index, alpha, result);
    }

    ArBBKernels5D::ArBBKernels5D() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;
    }

    ArBBKernels5D::~ArBBKernels5D() {
    }

    double ArBBKernels5D::multTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_source;

        if (isMultTransfirst && isMultfirst) {
          ArBB_Level_5D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_5D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_5D.write_only_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_5D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1], ptrLevel[(i * dims) + 2], ptrLevel[(i * dims) + 3], ptrLevel[(i * dims) + 4]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1], ptrIndex[(i * dims) + 2], ptrIndex[(i * dims) + 3], ptrIndex[(i * dims) + 4]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_5D = arbb::dense< vecElem64 >(sourceSize);

          arbb::range< vecElem64 > Data_range = ArBB_Data_5D.read_write_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1], ptrData[(i * dims) + 2], ptrData[(i * dims) + 3], ptrData[(i * dims) + 4]};

            Data_range[i] = curData;
          }

          isMultfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f64>))(ArBB_Data_5D, ArBB_Level_5D, ArBB_Index_5D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels5D::multArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_alpha;

        if (isMultTransfirst && isMultfirst) {
          ArBB_Level_5D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_5D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_5D.read_write_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_5D.read_write_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1], ptrLevel[(i * dims) + 2], ptrLevel[(i * dims) + 3], ptrLevel[(i * dims) + 4]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1], ptrIndex[(i * dims) + 2], ptrIndex[(i * dims) + 3], ptrIndex[(i * dims) + 4]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_5D = arbb::dense< vecElem64 >(result_size);

          arbb::range< vecElem64 > Data_range = ArBB_Data_5D.read_write_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1], ptrData[(i * dims) + 2], ptrData[(i * dims) + 3], ptrData[(i * dims) + 4]};

            Data_range[i] = curData;
          }

          isMultTransfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f64>))(ArBB_Data_5D, ArBB_Level_5D, ArBB_Index_5D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels5D::multTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_source;

        if (isMultTransSPfirst && isMultSPfirst) {
          ArBB_LevelSP_5D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_5D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_5D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_5D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1], ptrLevel[(i * dims) + 2], ptrLevel[(i * dims) + 3], ptrLevel[(i * dims) + 4]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1], ptrIndex[(i * dims) + 2], ptrIndex[(i * dims) + 3], ptrIndex[(i * dims) + 4]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_5D = arbb::dense< vecElem32 >(sourceSize);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_5D.write_only_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1], ptrData[(i * dims) + 2], ptrData[(i * dims) + 3], ptrData[(i * dims) + 4]};

            DataSP_range[i] = curData;
          }

          isMultSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f32>))(ArBB_DataSP_5D, ArBB_LevelSP_5D, ArBB_IndexSP_5D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels5D::multSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      try {
        arbb::dense<arbb::f32, 1> ArBB_result;
        arbb::dense<arbb::f32, 1> ArBB_alpha;

        if (isMultTransSPfirst && isMultSPfirst) {
          ArBB_LevelSP_5D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_5D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_5D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_5D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1], ptrLevel[(i * dims) + 2], ptrLevel[(i * dims) + 3], ptrLevel[(i * dims) + 4]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1], ptrIndex[(i * dims) + 2], ptrIndex[(i * dims) + 3], ptrIndex[(i * dims) + 4]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_5D = arbb::dense< vecElem32 >(result_size);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_5D.write_only_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1], ptrData[(i * dims) + 2], ptrData[(i * dims) + 3], ptrData[(i * dims) + 4]};

            DataSP_range[i] = curData;
          }

          isMultTransSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f32>))(ArBB_DataSP_5D, ArBB_LevelSP_5D, ArBB_IndexSP_5D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    void ArBBKernels5D::resetKernels() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;
    }

  }

}
