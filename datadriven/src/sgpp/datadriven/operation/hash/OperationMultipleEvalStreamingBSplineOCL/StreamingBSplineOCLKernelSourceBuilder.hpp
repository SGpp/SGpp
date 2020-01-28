// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <fstream>
#include <string>

#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace datadriven {

namespace streamingBSplineOCL {
template <typename T>
struct getType {};

template <>
struct getType<float> {
  static std::string asString() { return "float"; }
  static std::string constSuffix() { return "f"; }
  static std::string intAsString() { return "uint"; }
};

template <>
struct getType<double> {
  static std::string asString() { return "double"; }
  static std::string constSuffix() { return ""; }
  static std::string intAsString() { return "ulong"; }
};
}  // namespace streamingBSplineOCL

template <typename real_type>
class StreamingBSplineOCLKernelSourceBuilder {
 private:
  size_t degree;
  std::shared_ptr<base::OCLOperationConfiguration> parameters;

 public:
  StreamingBSplineOCLKernelSourceBuilder(
      size_t degree, std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : degree(degree), parameters(parameters) {}

  std::string generateSourceMult(size_t dims) {
    if ((*parameters)["REUSE_SOURCE"].getBool()) {
      std::stringstream streamProgramSrc;
      std::ifstream file;
      file.open("multKernelBSpline_tmp.cl");

      if (file.is_open()) {
        std::string line;

        while (getline(file, line)) {
          streamProgramSrc << line << std::endl;
        }

        file.close();
      } else {
        throw base::operation_exception("OCL error: file to reuse not found\n");
      }

      return streamProgramSrc.str();
    }

    size_t localWorkgroupSize = (*parameters)["LOCAL_SIZE"].getUInt();
    bool useLocalMemory = (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool();

    size_t dataBlockSize = (*parameters)["KERNEL_DATA_BLOCK_SIZE"].getUInt();

    std::string streamProgramSrc;
    const std::string dp1h = getDegreePlusOneHalvedString();
    const std::string evalFormula = getBSplineEvalFormula();

    if (streamingBSplineOCL::getType<real_type>::asString() == "double") {
      streamProgramSrc += R"(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
)";
    }

    streamProgramSrc +=
        R"(
__kernel
__attribute__((reqd_work_group_size({localWorkgroupSize}, 1, 1)))
void multOCL(__global const {float}* ptrLevel,
__global const {float}* ptrIndex,
__global const {float}* ptrData,
__global const {float}* ptrAlpha,
__global       {float}* ptrResult,
uint resultSize,
uint start_grid,
uint end_grid) {
  int globalIdx = get_global_id(0);
  int localIdx = get_local_id(0);
  {float} tmp;
)";

    for (size_t i = 0; i < degree; i++) {
      streamProgramSrc += replace(R"(
  {float} c{i};
)",
                                  "{i}", i);
    }

    if (useLocalMemory) {
      streamProgramSrc +=
          R"(
  __local {float} locLevel[{dimsTimesLocalWorkgroupSize}];
  __local {float} locIndex[{dimsTimesLocalWorkgroupSize}];
  __local {float} locAlpha[{localWorkgroupSize}];
)";
    }

    for (size_t i = 0; i < dataBlockSize; i++) {
      streamProgramSrc += replace(R"(
  {float} curSupport_{i};
  {float} myResult_{i} = 0.0;
)",
                                  "{i}", i);
    }

    streamProgramSrc += R"(
  // create registers for the data)";

    // TODO(pfandedd): might lead to bad access pattern, as each 1D array is accessed with a stride
    // of dim
    // |***|***|***|*** -> better: ||||************** and then ****||||************
    for (size_t i = 0; i < dataBlockSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        streamProgramSrc += replace(
            R"(
  {float} data_{i}_{d} = ptrData[{i} + ({dataBlockSize} * globalIdx) +
                                 (resultSize * {d})];
  {float} x_{i}_{d};
)",
            "{d}", d);
      }

      replaceInPlace(streamProgramSrc, "{i}", i);
      streamProgramSrc += "\n";
    }

    streamProgramSrc += "\n";

    if (useLocalMemory) {
      streamProgramSrc +=
          R"(
  // iterate over all grid points (fast ones, with cache)
  uint chunkSizeGrid = end_grid - start_grid;
  uint fastChunkSizeGrid = (chunkSizeGrid / {localWorkgroupSize}) *
                           {localWorkgroupSize};
  for (int j = start_grid; j < start_grid + fastChunkSizeGrid;
       j += {localWorkgroupSize}) {
)";

      for (size_t d = 0; d < dims; d++) {
        streamProgramSrc += replace(
            R"(
    locLevel[(localIdx*{dims})+{d}] = ptrLevel[((j+localIdx)*{dims})+{d}];
    locIndex[(localIdx*{dims})+{d}] = ptrIndex[((j+localIdx)*{dims})+{d}];)",
            "{d}", d);
      }

      streamProgramSrc +=
          R"(
    locAlpha[localIdx] = ptrAlpha[j+localIdx];
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int k = 0; k < {localWorkgroupSize}; k++) {
)";

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
      curSupport_{i} = locAlpha[k];)",
                                    "{i}", i);

        for (size_t d = 0; d < dims; d++) {
          streamProgramSrc += replace(replace(replace(
                                                  R"(
    x_{i}_{d} = data_{i}_{d} * locLevel[(k*{dims})+{d}] - 
                locIndex[(k*{dims})+{d}] + {dp1h};)",
                                                  "{i}", i),
                                              "{d}", d),
                                      "{dp1h}", dp1h);
        }
      }

      for (size_t d = 0; d < dims; d++) {
        /*streamProgramSrc += R"(
         if ((locLevel[(k*{dims})+{d}]) == 2.0{f}) {
         } else if ((locIndex[(k*{dims})+{d}]) == 1.0{f}) {
         )";*/

        for (size_t i = 0; i < dataBlockSize; i++) {
          streamProgramSrc += replace(replace(evalFormula, "{x}", "x_{i}_{d}"), "{i}", i);

          /*streamProgramSrc += replace(R"(
           curSupport_{i} *= max(2.0{f} - ((locLevel[(k*{dims})+{d}]) *
           (data_{i}_{d})), 0.0{f});)",
           "{i}", i);*/
        }

        /*streamProgramSrc += R"(
         } else if ((locIndex[(k*{dims})+{d}]) ==
         ((locLevel[(k*{dims})+{d}]) - 1.0{f})) {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(((locLevel[(k*{dims})+{d}]) * (data_{i}_{d})) -
         (locIndex[(k*{dims})+{d}]) + 1.0{f}, 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         } else {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(
         1.0{f} - fabs(((locLevel[(k*{dims})+{d}]) * (data_{i}_{d})) -
         (locIndex[(k*{dims})+{d}])), 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         }
         )";*/

        replaceInPlace(streamProgramSrc, "{d}", d);
      }

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
      myResult_{i} += curSupport_{i};)",
                                    "{i}", i);
      }

      streamProgramSrc +=
          R"(
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  // iterate over all grid points (slow ones, without cache)
  for (int m = start_grid + fastChunkSizeGrid; m < end_grid; m++) {
)";

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
    curSupport_{i} = ptrAlpha[m];)",
                                    "{i}", i);

        for (size_t d = 0; d < dims; d++) {
          streamProgramSrc += replace(replace(replace(
                                                  R"(
    x_{i}_{d} = data_{i}_{d} * ptrLevel[(m*{dims})+{d}] - 
                ptrIndex[(m*{dims})+{d}] + {dp1h};)",
                                                  "{i}", i),
                                              "{d}", d),
                                      "{dp1h}", dp1h);
        }
      }

      for (size_t d = 0; d < dims; d++) {
        /*streamProgramSrc += R"(
         if ((ptrLevel[(m*{dims})+{d}]) == 2.0{f}) {
         } else if ((ptrIndex[(m*{dims})+{d}]) == 1.0{f}) {
         )";*/

        for (size_t i = 0; i < dataBlockSize; i++) {
          streamProgramSrc += replace(replace(evalFormula, "{x}", "x_{i}_{d}"), "{i}", i);

          /*streamProgramSrc += replace(R"(
           curSupport_{i} *= max(2.0{f} - ((ptrLevel[(m*{dims})+{d}]) *
           (data_{i}_{d})), 0.0{f});)",
           "{i}", i);*/
        }

        /*streamProgramSrc += R"(
         } else if ((ptrIndex[(m*{dims})+{d}]) ==
         ((ptrLevel[(m*{dims})+{d}]) - 1.0{f})) {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(((ptrLevel[(m*{dims})+{d}]) * (data_{i}_{d})) -
         (ptrIndex[(m*{dims})+{d}]) + 1.0{f}, 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         } else {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(
         1.0{f} - fabs(((ptrLevel[(m*{dims})+{d}]) * (data_{i}_{d})) -
         (ptrIndex[(m*{dims})+{d}])), 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         }
         )";*/

        replaceInPlace(streamProgramSrc, "{d}", d);
      }

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
    myResult_{i} += curSupport_{i};)",
                                    "{i}", i);
      }

      streamProgramSrc += R"(
  }
)";
    } else {
      streamProgramSrc +=
          R"(
  // iterate over all grid points (slow ones, without cache)
  for (int m = start_grid; m < end_grid; m++) {
)";

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
    curSupport_{i} = ptrAlpha[m];)",
                                    "{i}", i);

        for (size_t d = 0; d < dims; d++) {
          streamProgramSrc += replace(replace(replace(
                                                  R"(
    x_{i}_{d} = data_{i}_{d} * ptrLevel[(m*{dims})+{d}] - 
                ptrIndex[(m*{dims})+{d}] + {dp1h};)",
                                                  "{i}", i),
                                              "{d}", d),
                                      "{dp1h}", dp1h);
        }
      }

      for (size_t d = 0; d < dims; d++) {
        /*streamProgramSrc += R"(
         if ((ptrLevel[(m*{dims})+{d}]) == 2.0{f}) {
         } else if ((ptrIndex[(m*{dims})+{d}]) == 1.0{f}) {
         )";*/

        for (size_t i = 0; i < dataBlockSize; i++) {
          streamProgramSrc += replace(replace(evalFormula, "{x}", "x_{i}_{d}"), "{i}", i);

          /*streamProgramSrc += replace(R"(
           curSupport_{i} *= max(2.0{f} - ((ptrLevel[(m*{dims})+{d}]) *
           (data_{i}_{d})), 0.0{f});)",
           "{i}", i);*/
        }

        /*streamProgramSrc += R"(
         } else if ((ptrIndex[(m*{dims})+{d}]) ==
         ((ptrLevel[(m*{dims})+{d}]) - 1.0{f})) {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(((ptrLevel[(m*{dims})+{d}]) * (data_{i}_{d})) -
         (ptrIndex[(m*{dims})+{d}]) + 1.0{f}, 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         } else {
         )";

         for (size_t i = 0; i < dataBlockSize; i++) {
         streamProgramSrc += replace(R"(
         curSupport_{i} *= max(
         1.0{f} - fabs(((ptrLevel[(m*{dims})+{d}]) * (data_{i}_{d})) -
         (ptrIndex[(m*{dims})+{d}])), 0.0{f});)",
         "{i}", i);
         }

         streamProgramSrc += R"(
         }
         )";*/
        replaceInPlace(streamProgramSrc, "{d}", d);
      }

      for (size_t i = 0; i < dataBlockSize; i++) {
        streamProgramSrc += replace(R"(
    myResult_{i} += curSupport_{i};)",
                                    "{i}", i);
      }

      streamProgramSrc += R"(
  }
)";
    }

    streamProgramSrc += "\n";

    for (size_t i = 0; i < dataBlockSize; i++) {
      //      streamProgramSrc << "         printf(\"myResult_" << i << ": %lf\\n\", myResult_" << i
      //      << ");\n";
      //      streamProgramSrc << "         printf(\"writing to index: %i\\n\", (" << dataBlockSize
      //      << " * globalIdx) + " << i
      //          << ");\n";
      streamProgramSrc += replace(R"(
  ptrResult[({dataBlockSize} * globalIdx) + {i}] = myResult_{i};)",
                                  "{i}", i);
    }

    streamProgramSrc += R"(
}
)";

    replaceInPlace(streamProgramSrc, "{float}",
                   streamingBSplineOCL::getType<real_type>::asString());
    replaceInPlace(streamProgramSrc, "{f}", streamingBSplineOCL::getType<real_type>::constSuffix());
    replaceInPlace(streamProgramSrc, "{dataBlockSize}", dataBlockSize);
    replaceInPlace(streamProgramSrc, "{dims}", dims);
    replaceInPlace(streamProgramSrc, "{localWorkgroupSize}", localWorkgroupSize);
    replaceInPlace(streamProgramSrc, "{dimsTimesLocalWorkgroupSize}", dims * localWorkgroupSize);

    // update file with kernel (for debugging)
    std::ofstream multFile;
    multFile.open("multKernelBSpline_tmp.cl");
    multFile << streamProgramSrc;
    multFile.close();

    return streamProgramSrc;
  }

  std::string generateSourceMultTrans(size_t dims) {
    if ((*parameters)["REUSE_SOURCE"].getBool()) {
      std::stringstream streamProgramSrc;
      std::ifstream file;
      file.open("multTransKernelBSpline_tmp.cl");

      if (file.is_open()) {
        std::string line;

        while (getline(file, line)) {
          streamProgramSrc << line << std::endl;
        }

        file.close();
      } else {
        throw base::operation_exception("OCL error: file to reuse not found\n");
      }

      return streamProgramSrc.str();
    }

    size_t localWorkgroupSize = (*parameters)["LOCAL_SIZE"].getUInt();
    size_t transDataBlockSize = (*parameters)["KERNEL_TRANS_DATA_BLOCK_SIZE"].getUInt();

    std::string streamProgramSrc;
    const std::string dp1h = getDegreePlusOneHalvedString();
    const std::string evalFormula = getBSplineEvalFormula();

    if (streamingBSplineOCL::getType<real_type>::asString() == "double") {
      streamProgramSrc += R"(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
)";
    }

    streamProgramSrc +=
        R"(
__kernel
__attribute__((reqd_work_group_size({localWorkgroupSize}, 1, 1)))
void multTransOCL(__global const {float}* ptrLevel,
__global const {float}* ptrIndex,
__global const {float}* ptrData,
__global const {float}* ptrSource,
__global       {float}* ptrResult,
uint sourceSize,
uint start_data,
uint end_data) {
  int globalIdx = get_global_id(0);
  int groupIdx = get_group_id(0);
  int localIdx = get_local_id(0);

  __local double resultsTemp[{localWorkgroupSize}];

  {float} tmp;
  {float} myResult = 0.0;
)";

    for (size_t i = 0; i < degree; i++) {
      streamProgramSrc += replace(R"(
  {float} c{i};
)",
                                  "{i}", i);
    }

    for (size_t d = 0; d < dims; d++) {
      streamProgramSrc += replace(
          R"(
  {float} level_{d} = ptrLevel[(groupIdx*{dims})+{d}];
  {float} index_{d} = ptrIndex[(groupIdx*{dims})+{d}];)",
          "{d}", d);
    }

    streamProgramSrc += replace(R"(

  for (int k = start_data + localIdx; k < end_data; k += {upper_bound}) {
)",
                                "{upper_bound}", transDataBlockSize * localWorkgroupSize);

    for (size_t i = 0; i < transDataBlockSize; i++) {
      streamProgramSrc += replace(replace(R"(
    {float} curSupport_{i} = ptrSource[k + {iTimesLWGS}];)",
                                          "{i}", i),
                                  "{iTimesLWGS}", localWorkgroupSize * i);
    }

    // streamProgramSrc += "\n";

    for (size_t i = 0; i < transDataBlockSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        streamProgramSrc += replace(replace(
                                        R"(
    {float} x_{i}_{d} = ptrData[({d}*sourceSize)+k + {iTimesLWGS}] * level_{d}
                        - index_{d} + {dp1h};)",
                                        "{d}", d),
                                    "{dp1h}", dp1h);
      }

      streamProgramSrc =
          replace(replace(streamProgramSrc, "{i}", i), "{iTimesLWGS}", localWorkgroupSize * i);
      streamProgramSrc += "\n";
    }

    for (size_t d = 0; d < dims; d++) {
      /*streamProgramSrc += R"(
       if ((level_{d}) == 2.0{f}) {
       } else if ((index_{d}) = 1.0{f}) {
       )";*/

      for (size_t i = 0; i < transDataBlockSize; i++) {
        streamProgramSrc += replace(replace(evalFormula, "{x}", "x_{i}_{d}"), "{i}", i);

        /*streamProgramSrc += replace(replace(R"(
         curSupport_{i} *= max(
         2.0{f} - ((level_{d}) * (ptrData[({d}*sourceSize)+k +
         {iTimesLWGS}])), 0.0{f});)",
         "{i}", i), "{iTimesLWGS}", localWorkgroupSize * i);*/
      }

      /*streamProgramSrc += R"(
       } else if ((index_{d}) == ((level_{d}) - 1.0{f})) {
       )";

       for (size_t i = 0; i < transDataBlockSize; i++) {
       streamProgramSrc += replace(replace(
       R"(
       curSupport_{i} *= max(
       ((level_{d}) * (ptrData[({d}*sourceSize)+k + {iTimesLWGS}])) -
       (index_{d}) + 1.0{f}, 0.0{f});)",
       "{i}", i), "{iTimesLWGS}", localWorkgroupSize * i);
       }

       streamProgramSrc += R"(
       } else {
       )";

       for (size_t i = 0; i < transDataBlockSize; i++) {
       streamProgramSrc += replace(replace(R"(
       curSupport_{i} *= max(
       1.0{f} - fabs(((level_{d}) *
       (ptrData[({d}*sourceSize)+k + {iTimesLWGS}])) -
       (index_{d})), 0.0{f});)",
       "{i}", i), "{iTimesLWGS}", localWorkgroupSize * i);
       }

       streamProgramSrc +=  R"(
       }
       )";*/

      replaceInPlace(streamProgramSrc, "{d}", d);
    }

    for (size_t i = 0; i < transDataBlockSize; i++) {
      streamProgramSrc += replace(R"(
    myResult += curSupport_{i};)",
                                  "{i}", i);
    }

    streamProgramSrc +=
        R"(
  }

  resultsTemp[localIdx] = myResult;
  barrier(CLK_LOCAL_MEM_FENCE);

  if (localIdx == 0) {
    double overallResult = 0.0;
    for (int i = 0; i < {localWorkgroupSize}; i++) {
      overallResult += resultsTemp[i];
    }
    ptrResult[groupIdx] = overallResult;
  }
}
)";

    replaceInPlace(streamProgramSrc, "{float}",
                   streamingBSplineOCL::getType<real_type>::asString());
    replaceInPlace(streamProgramSrc, "{f}", streamingBSplineOCL::getType<real_type>::constSuffix());
    replaceInPlace(streamProgramSrc, "{dims}", dims);
    replaceInPlace(streamProgramSrc, "{localWorkgroupSize}", localWorkgroupSize);

    // update file with kernel (for debugging)
    std::ofstream multTransFile;
    multTransFile.open("multTransKernelBSpline_tmp.cl");
    multTransFile << streamProgramSrc;
    multTransFile.close();

    return streamProgramSrc;
  }

 protected:
  std::string getDegreePlusOneHalvedString() const {
    return std::to_string((degree + 1) / 2) + ((degree % 2 == 0) ? ".5{f}" : ".0{f}");
  }

  std::string getBSplineEvalFormula() const {
    switch (degree) {
      case 1:
        return R"(
    if (({x} < 0.0{f}) || ({x} >= 2.0{f})) {
      curSupport_{i} = 0.0{f};
    } else if ({x} < 1.0{f}) {
      curSupport_{i} *= {x};
    } else {
      curSupport_{i} *= 2.0{f} - {x};
    }
)";
        break;
      case 3:
        return R"(
    if (({x} < 0.0{f}) || ({x} >= 4.0{f})) {
      curSupport_{i} = 0.0{f};
    } else {
      if ({x} < 1.0{f}) {
        tmp = 1.0{f} / 6.0{f};
        c2 = 0.0{f};
        c1 = 0.0{f};
        c0 = 0.0{f};
      } else if ({x} < 2.0{f}) {
        tmp = -0.5{f};
        c2 = 2.0{f};
        c1 = -2.0{f};
        c0 = 2.0{f} / 3.0{f};
      } else if ({x} < 3.0{f}) {
        tmp = 0.5{f};
        c2 = -4.0{f};
        c1 = 10.0{f};
        c0 = -22.0{f} / 3.0{f};
      } else {
        tmp = -1.0{f} / 6.0{f};
        c2 = 2.0{f};
        c1 = -8.0{f};
        c0 = 32.0{f} / 3.0{f};
      }

      tmp = c2 + tmp * {x};
      tmp = c1 + tmp * {x};
      tmp = c0 + tmp * {x};
      curSupport_{i} *= tmp;
    }
)";
        /*return R"(
         if (({x} < 0.0{f}) || ({x} >= 4.0{f})) {
         curSupport_{i} = 0.0{f};
         } else if ({x} < 1.0{f}) {
         tmp = 1.0{f} / 6.0{f};
         tmp *= {x};
         tmp *= {x};
         tmp *= {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 2.0{f}) {
         tmp = -0.5{f};
         tmp = 2.0{f} + tmp * {x};
         tmp = -2.0{f} + tmp * {x};
         tmp = (2.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 3.0{f}) {
         tmp = 0.5{f};
         tmp = -4.0{f} + tmp * {x};
         tmp = 10.0{f} + tmp * {x};
         tmp = -(22.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         } else {
         tmp = -(1.0{f} / 6.0{f});
         tmp = 2.0{f} + tmp * {x};
         tmp = -8.0{f} + tmp * {x};
         tmp = (32.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         }
         )";*/
        // this takes longer despite fewer if-clauses on average
        // (assuming the cases are reached with the same probability)
        /*return R"(
         if (({x} < 0.0{f}) || ({x} >= 4.0{f})) {
         curSupport_{i} = 0.0{f};
         } else if ({x} < 2.0{f}) {
         if ({x} < 1.0{f}) {
         tmp = 1.0{f} / 6.0{f};
         tmp *= {x};
         tmp *= {x};
         tmp *= {x};
         curSupport_{i} *= tmp;
         } else {
         tmp = -0.5{f};
         tmp = 2.0{f} + tmp * {x};
         tmp = -2.0{f} + tmp * {x};
         tmp = (2.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         }
         } else {
         if ({x} < 3.0{f}) {
         tmp = 0.5{f};
         tmp = -4.0{f} + tmp * {x};
         tmp = 10.0{f} + tmp * {x};
         tmp = -(22.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         } else {
         tmp = -(1.0{f} / 6.0{f});
         tmp = 2.0{f} + tmp * {x};
         tmp = -8.0{f} + tmp * {x};
         tmp = (32.0{f} / 3.0{f}) + tmp * {x};
         curSupport_{i} *= tmp;
         }
         }
         )";*/
        break;
      case 5:
        return R"(
    if (({x} < 0.0{f}) || ({x} >= 6.0{f})) {
      curSupport_{i} = 0.0{f};
    } else {
      if ({x} < 1.0{f}) {
        tmp = 1.0{f} / 120.0{f};
        c4 = 0.0{f};
        c3 = 0.0{f};
        c2 = 0.0{f};
        c1 = 0.0{f};
        c0 = 0.0{f};
      } else if ({x} < 2.0{f}) {
        tmp = -1.0{f} / 24.0{f};
        c4 = 0.25{f};
        c3 = -0.5{f};
        c2 = 0.5{f};
        c1 = -0.25{f};
        c0 = 0.05{f};
      } else if ({x} < 3.0{f}) {
        tmp = 1.0{f} / 12.0{f};
        c4 = -1.0{f};
        c3 = 4.5{f};
        c2 = -9.5{f};
        c1 = 9.75{f};
        c0 = -3.95{f};
      } else if ({x} < 4.0{f}) {
        tmp = -1.0{f} / 12.0{f};
        c4 = 1.5{f};
        c3 = -10.5{f};
        c2 = 35.5{f};
        c1 = -57.75{f};
        c0 = 36.55{f};
      } else if ({x} < 5.0{f}) {
        tmp = 1.0{f} / 24.0{f};
        c4 = -1.0{f};
        c3 = 9.5{f};
        c2 = -44.5{f};
        c1 = 102.25{f};
        c0 = -91.45{f};
      } else {
        tmp = -1.0{f} / 120.0{f};
        c4 = 0.25{f};
        c3 = -3.0{f};
        c2 = 18.0{f};
        c1 = -54.0{f};
        c0 = 64.8{f};
      }

      tmp = c4 + tmp * {x};
      tmp = c3 + tmp * {x};
      tmp = c2 + tmp * {x};
      tmp = c1 + tmp * {x};
      tmp = c0 + tmp * {x};
      curSupport_{i} *= tmp;
    }
)";
        /*return R"(
         if (({x} < 0.0{f}) || ({x} >= 6.0{f})) {
         curSupport_{i} = 0.0{f};
         } else if ({x} < 1.0{f}) {
         tmp = 1.0{f} / 120.0{f};
         tmp *= {x};
         tmp *= {x};
         tmp *= {x};
         tmp *= {x};
         tmp *= {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 2.0{f}) {
         tmp = -1.0{f} / 24.0{f};
         tmp = 0.25{f} + tmp * {x};
         tmp = -0.5{f} + tmp * {x};
         tmp = 0.5{f} + tmp * {x};
         tmp = -0.25{f} + tmp * {x};
         tmp = 0.05{f} + tmp * {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 3.0{f}) {
         tmp = 1.0{f} / 12.0{f};
         tmp = -1.0{f} + tmp * {x};
         tmp = 4.5{f} + tmp * {x};
         tmp = -9.5{f} + tmp * {x};
         tmp = 9.75{f} + tmp * {x};
         tmp = -3.95{f} + tmp * {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 4.0{f}) {
         tmp = -1.0{f} / 12.0{f};
         tmp = 1.5{f} + tmp * {x};
         tmp = -10.5{f} + tmp * {x};
         tmp = 35.5{f} + tmp * {x};
         tmp = -57.75{f} + tmp * {x};
         tmp = 36.55{f} + tmp * {x};
         curSupport_{i} *= tmp;
         } else if ({x} < 5.0{f}) {
         tmp = 1.0{f} / 24.0{f};
         tmp = -1.0{f} + tmp * {x};
         tmp = 9.5{f} + tmp * {x};
         tmp = -44.5{f} + tmp * {x};
         tmp = 102.25{f} + tmp * {x};
         tmp = -91.45{f} + tmp * {x};
         curSupport_{i} *= tmp;
         } else {
         tmp = -1.0{f} / 120.0{f};
         tmp = 0.25{f} + tmp * {x};
         tmp = -3.0{f} + tmp * {x};
         tmp = 18.0{f} + tmp * {x};
         tmp = -54.0{f} + tmp * {x};
         tmp = 64.8{f} + tmp * {x};
         curSupport_{i} *= tmp;
         }
         )";*/
        break;
      default:
        throw base::operation_exception("degree not supported.");
    }
  }

  static void replaceInPlace(std::string& str, const std::string& searchStr,
                             const std::string& replaceStr) {
    size_t pos = 0;

    while ((pos = str.find(searchStr, pos)) != std::string::npos) {
      str.replace(pos, searchStr.length(), replaceStr);
      pos += replaceStr.length();
    }
  }

  static void replaceInPlace(std::string& str, const std::string& searchStr, size_t number) {
    replaceInPlace(str, searchStr, std::to_string(number));
  }

  static std::string replace(const std::string& str, const std::string& searchStr,
                             const std::string& replaceStr) {
    std::string result = str;
    replaceInPlace(result, searchStr, replaceStr);
    return result;
  }

  static std::string replace(const std::string& str, const std::string& searchStr,
                             size_t replaceNumber) {
    return replace(str, searchStr, std::to_string(replaceNumber));
  }
};
}  // namespace datadriven
}  // namespace sgpp
