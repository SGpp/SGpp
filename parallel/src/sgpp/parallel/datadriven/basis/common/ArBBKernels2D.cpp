/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/datadriven/basis/common/ArBBKernels2D.hpp>

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <string>

#include <arbb.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace parallel {

    size_t global_dims_2D;

    typedef arbb::array<arbb::f64, 2> vecElem64;
    typedef arbb::uncaptured< arbb::array<double, 2> >::type uc_vecElem64;
    typedef arbb::array<arbb::f32, 2> vecElem32;
    typedef arbb::uncaptured< arbb::array<float, 2> >::type uc_vecElem32;

    arbb::dense< vecElem32, 1> ArBB_DataSP_2D;
    arbb::dense< vecElem32, 1> ArBB_LevelSP_2D;
    arbb::dense< vecElem32, 1> ArBB_IndexSP_2D;

    arbb::dense< vecElem64, 1> ArBB_Data_2D;
    arbb::dense< vecElem64, 1> ArBB_Level_2D;
    arbb::dense< vecElem64, 1> ArBB_Index_2D;

    template <typename fp_Type>
    void arbb_multTrans(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPointTrans(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::array<fp_Type, 2>& Level, const arbb::array<fp_Type, 2>& Index, const arbb::dense<fp_Type>& source, fp_Type& result) {
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
    void arbb_multTrans_mod(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPointTrans(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::array<fp_Type, 2>& Level, const arbb::array<fp_Type, 2>& Index, const arbb::dense<fp_Type>& source, fp_Type& result) {
          result = 0.0;

          _for (arbb::usize i = 0, i < Data.num_cols(), i++) {
            arbb::array<fp_Type, 2> curData = Data[i];
            fp_Type tmp = source[i];

            for (size_t l = 0; l < 2; l++) {
              _if (Level[l] == 2.0) {
              }
              _else_if (Index[l] == 1.0) {
                tmp *= arbb::max(2.0 - (Level[l] * curData[l]), 0.0);
              }
              _else_if (Index[l] == (Level[l] - 1.0)) {
                tmp *= arbb::max((Level[l] * curData[l]) - Index[l] + 1.0, 0.0);
              }
              _else {
                tmp *= arbb::max((1.0 - arbb::abs(((Level[l] * curData[l]) - Index[l]))), 0.0);
              } _end_if;
            }

            result += tmp;
          }
          _end_for;
        }
      };

      arbb::map(local::evalGridPointTrans)(Data, Level, Index, source, result);
    }

    template <typename fp_Type>
    void arbb_mult(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPoint(const arbb::array<fp_Type, 2>& DataPoint, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& alpha, fp_Type& result) {
          result = 0.0;

          _for (arbb::usize j = 0, j < Level.num_cols(), j++) {
            result += (alpha[j] * arbb::mul_reduce(arbb::max((1.0 - arbb::abs(((Level[j] * DataPoint) - Index[j]))), 0.0)));
          }
          _end_for;
        }
      };

      arbb::map(local::evalGridPoint)(Data, Level, Index, alpha, result);
    }

    template <typename fp_Type>
    void arbb_mult_mod(const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Data, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result) {
      struct local {
        static void evalGridPoint(const arbb::array<fp_Type, 2>& DataPoint, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Level, const arbb::dense< arbb::array<fp_Type, 2>, 1 >& Index, const arbb::dense<fp_Type>& alpha, fp_Type& result) {
          result = 0.0;

          _for (arbb::usize j = 0, j < Level.num_cols(), j++) {
            fp_Type tmp = alpha[j];
            arbb::array<fp_Type, 2> curLevel = Level[j];
            arbb::array<fp_Type, 2> curIndex = Index[j];

            for (size_t l = 0; l < 2; l++) {
              _if (curLevel[l] == 2.0) {
              }
              _else_if (curIndex[l] == 1.0) {
                tmp *= arbb::max(2.0 - (curLevel[l] * DataPoint[l]), 0.0);
              }
              _else_if (curIndex[l] == (curLevel[l] - 1.0)) {
                tmp *= arbb::max((curLevel[l] * DataPoint[l]) - curIndex[l] + 1.0, 0.0);
              }
              _else {
                tmp *= arbb::max((1.0 - arbb::abs(((curLevel[l] * DataPoint[l]) - curIndex[l]))), 0.0);
              } _end_if;
            }

            result += tmp;
          }
          _end_for;
        }
      };

      arbb::map(local::evalGridPoint)(Data, Level, Index, alpha, result);
    }

    ArBBKernels2D::ArBBKernels2D() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;

      isMultModTransSPfirst = true;
      isMultModSPfirst = true;

      isMultModTransfirst = true;
      isMultModfirst = true;
    }

    ArBBKernels2D::~ArBBKernels2D() {
    }

    double ArBBKernels2D::multTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_source;

        if (isMultTransfirst && isMultfirst) {
          ArBB_Level_2D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_2D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_2D.write_only_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_2D = arbb::dense< vecElem64 >(sourceSize);

          arbb::range< vecElem64 > Data_range = ArBB_Data_2D.read_write_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            Data_range[i] = curData;
          }

          isMultfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f64>))(ArBB_Data_2D, ArBB_Level_2D, ArBB_Index_2D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_alpha;

        if (isMultTransfirst && isMultfirst) {
          ArBB_Level_2D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_2D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_2D.read_write_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_2D.read_write_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_2D = arbb::dense< vecElem64 >(result_size);

          arbb::range< vecElem64 > Data_range = ArBB_Data_2D.read_write_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            Data_range[i] = curData;
          }

          isMultTransfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f64>))(ArBB_Data_2D, ArBB_Level_2D, ArBB_Index_2D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_source;

        if (isMultTransSPfirst && isMultSPfirst) {
          ArBB_LevelSP_2D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_2D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_2D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_2D = arbb::dense< vecElem32 >(sourceSize);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_2D.write_only_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            DataSP_range[i] = curData;
          }

          isMultSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f32>))(ArBB_DataSP_2D, ArBB_LevelSP_2D, ArBB_IndexSP_2D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f32, 1> ArBB_result;
        arbb::dense<arbb::f32, 1> ArBB_alpha;

        if (isMultTransSPfirst && isMultSPfirst) {
          ArBB_LevelSP_2D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_2D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_2D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_2D = arbb::dense< vecElem32 >(result_size);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_2D.write_only_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            DataSP_range[i] = curData;
          }

          isMultTransSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f32>))(ArBB_DataSP_2D, ArBB_LevelSP_2D, ArBB_IndexSP_2D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multModTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_source;

        if (isMultModTransfirst && isMultModfirst) {
          ArBB_Level_2D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_2D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_2D.write_only_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_2D = arbb::dense< vecElem64 >(sourceSize);

          arbb::range< vecElem64 > Data_range = ArBB_Data_2D.read_write_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            Data_range[i] = curData;
          }

          isMultModfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans_mod<arbb::f64>))(ArBB_Data_2D, ArBB_Level_2D, ArBB_Index_2D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multModArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_alpha;

        if (isMultModTransfirst && isMultModfirst) {
          ArBB_Level_2D = arbb::dense< vecElem64 >(storageSize);
          ArBB_Index_2D = arbb::dense< vecElem64 >(storageSize);

          arbb::range< vecElem64 > Level_range = ArBB_Level_2D.read_write_range();
          arbb::range< vecElem64 > Index_range = ArBB_Index_2D.read_write_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem64 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem64 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            Level_range[i] = curLevel;
            Index_range[i] = curIndex;
          }

          ArBB_Data_2D = arbb::dense< vecElem64 >(result_size);

          arbb::range< vecElem64 > Data_range = ArBB_Data_2D.read_write_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem64 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            Data_range[i] = curData;
          }

          isMultModTransfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult_mod<arbb::f64>))(ArBB_Data_2D, ArBB_Level_2D, ArBB_Index_2D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multModTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_source;

        if (isMultModTransSPfirst && isMultModSPfirst) {
          ArBB_LevelSP_2D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_2D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_2D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_2D = arbb::dense< vecElem32 >(sourceSize);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_2D.write_only_range();

          for (size_t i = 0; i < sourceSize; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            DataSP_range[i] = curData;
          }

          isMultModSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans_mod<arbb::f32>))(ArBB_DataSP_2D, ArBB_LevelSP_2D, ArBB_IndexSP_2D, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels2D::multModSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims_2D = dims;

      try {
        arbb::dense<arbb::f32, 1> ArBB_result;
        arbb::dense<arbb::f32, 1> ArBB_alpha;

        if (isMultModTransSPfirst && isMultModSPfirst) {
          ArBB_LevelSP_2D = arbb::dense< vecElem32 >(storageSize);
          ArBB_IndexSP_2D = arbb::dense< vecElem32 >(storageSize);

          arbb::range< vecElem32 > LevelSP_range = ArBB_LevelSP_2D.write_only_range();
          arbb::range< vecElem32 > IndexSP_range = ArBB_IndexSP_2D.write_only_range();

          for (size_t i = 0; i < storageSize; i++) {
            uc_vecElem32 curLevel = {ptrLevel[(i * dims) + 0], ptrLevel[(i * dims) + 1]};
            uc_vecElem32 curIndex = {ptrIndex[(i * dims) + 0], ptrIndex[(i * dims) + 1]};

            LevelSP_range[i] = curLevel;
            IndexSP_range[i] = curIndex;
          }

          ArBB_DataSP_2D = arbb::dense< vecElem32 >(result_size);

          arbb::range< vecElem32 > DataSP_range = ArBB_DataSP_2D.write_only_range();

          for (size_t i = 0; i < result_size; i++) {
            uc_vecElem32 curData = {ptrData[(i * dims) + 0], ptrData[(i * dims) + 1]};

            DataSP_range[i] = curData;
          }

          isMultModTransSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult_mod<arbb::f32>))(ArBB_DataSP_2D, ArBB_LevelSP_2D, ArBB_IndexSP_2D, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    void ArBBKernels2D::resetKernels() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;

      isMultModTransSPfirst = true;
      isMultModSPfirst = true;

      isMultModTransfirst = true;
      isMultModfirst = true;
    }

  }

}
