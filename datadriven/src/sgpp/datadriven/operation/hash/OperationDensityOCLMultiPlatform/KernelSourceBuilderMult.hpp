// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

#include <fstream>
#include <sstream>
#include <string>
namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// OpenCL source builder for density matrix vector multiplication
template <typename real_type>
class SourceBuilderMult : public base::KernelSourceBuilderBase<real_type> {
 private:
  /// OpenCL configuration containing the building flags
  json::Node &kernelConfiguration;
  /// Dimensions of grid
  size_t dims;
  /// Used workgroupsize for opencl kernel execution
  size_t localWorkgroupSize;
  /// Using local memory?
  bool useLocalMemory;
  size_t dataBlockSize;
  size_t transGridBlockSize;
  uint64_t maxDimUnroll;

  /// Use a cache for the 2^l values? Configuration parameter is USE_LEVEL_CACHE
  bool use_level_cache;
  /// Use a calculation scheme with less operations but more branching?
  bool use_less;
  /// Use ternary operator for branching? Configuration parameter is USE_LESS_OPERATIONS
  bool do_not_use_ternary;
  /// Avoid branching even at the cost of more operations? Configuration parameter is USE_IMPLICIT
  bool use_implicit_zero;
  /// Avoid using fmax? Configuration parameter is USE_FABS
  bool use_fabs_instead_of_fmax;
  /// Use preprocessed grid positions? Configuration parameter is PREPROCESSED_POSITIONS
  bool preprocess_positions;
  bool unroll_dim;

  /// Generate the opencl code to save the fixed gridpoint of a workitem to the local memory
  std::string save_from_global_to_private(size_t dimensions) {
    std::stringstream output;
    for (size_t block = 0; block < dataBlockSize; block++) {
      if (!preprocess_positions) {
        output << this->indent[0] << "__private int point_indices_block" << block << "["
               << dimensions << "];" << std::endl;
        output << this->indent[0] << "__private int point_level_block" << block << "[" << dimensions
               << "];" << std::endl;
        for (size_t i = 0; i < dimensions; i++) {
          output << this->indent[0] << "point_indices_block" << block << "[" << i
                 << "] = starting_points[(gridindex * " << dataBlockSize << " + " << block << ") * "
                 << dimensions << " * 2 + 2 * " << i << "];" << std::endl;
          output << this->indent[0] << "point_level_block" << block << "[" << i
                 << "] = starting_points[(gridindex * " << dataBlockSize << " + " << block << ") * "
                 << dimensions << " * 2 + 2 * " << i << " + 1];" << std::endl;
        }
      } else {
        output << this->indent[0] << "__private " << this->floatType() << " point_positions_block"
               << block << "[" << dimensions << "];" << std::endl;
        output << this->indent[0] << "__private " << this->floatType() << " point_hs_block" << block
               << "[" << dimensions << "];" << std::endl;
        output << this->indent[0] << "__private " << this->floatType() << " point_hinverses_block"
               << block << "[" << dimensions << "];" << std::endl;
        for (size_t i = 0; i < dimensions; i++) {
          // calculate inverses
          output << this->indent[1] << "point_hinverses_block" << block << "[" << i
                 << "] = hs_inverses[gridindex * " << dimensions << " + " << i << "];" << std::endl;
          // calculate hs
          output << this->indent[1] << "point_hs_block" << block << "[" << i
                 << "] = hs[gridindex * " << dimensions << " + " << i << "];" << std::endl;
          // calculate positions
          output << this->indent[1] << "point_positions_block" << block << "[" << i
                 << "] = positions[gridindex * " << dimensions << " + " << i << "];" << std::endl;
        }
      }
    }
    return output.str();
  }

  /// Generates the part of the opencl source code that calculates one entry of the density matrix
  std::string calculate_matrix_entry(size_t block, size_t dimensions) {
    std::stringstream output;
    // Use alias names for levels and indices
    std::string level_func1 =
        std::string("point_level_block") + std::to_string(block) + std::string("[dim]");
    std::string level_func2 = std::string("starting_points[i* ") + std::to_string(dimensions) +
                              std::string("*2+2*dim+1]");
    std::string index_func1 =
        std::string("point_indices_block") + std::to_string(block) + std::string("[dim]");
    std::string index_func2 =
        std::string("starting_points[i* ") + std::to_string(dimensions) + std::string("*2+2*dim]");
    // In case we use local memory we need to adjust the alias names
    if (useLocalMemory) {
      level_func2 =
          std::string("level_local[i* ") + std::to_string(dimensions) + std::string("+dim]");
      index_func2 =
          std::string("indices_local[i* ") + std::to_string(dimensions) + std::string("+dim]");
    }
    output << " zellenintegral = 1.0;" << std::endl;
    // In case we want replace the ternary operator we need to declare the counter variable
    if ((use_less && do_not_use_ternary) || (!use_implicit_zero && do_not_use_ternary))
      output << this->indent[2] << "int same_levels = " << dimensions << ";" << std::endl;
    // Loop over all dimensions
    output << this->indent[2] << "for(private int dim = 0;dim< " << dimensions << ";dim++) {"
           << std::endl;
    // In case we do not want to use that the entry is implicitly zero if we use the wrong order
    // we need to find the smallest level
    if (!use_implicit_zero) {
      // Initialise with one gridpoint
      output << this->indent[3] << "int index = " << index_func1 << ";" << std::endl;
      output << this->indent[3] << "int level = " << level_func1 << ";" << std::endl;
      output << this->indent[3] << "int index2 = " << index_func2 << ";" << std::endl;
      output << this->indent[3] << "int level2 = " << level_func2 << ";" << std::endl;
      // Check whether we need to use the other gridpoint as gridpoint 1
      output << this->indent[3] << "if (level > level2) {" << std::endl;
      output << this->indent[4] << "index = " << index_func2 << ";" << std::endl;
      output << this->indent[4] << "level = " << level_func2 << ";" << std::endl;
      output << this->indent[4] << "index2 = " << index_func1 << ";" << std::endl;
      output << this->indent[4] << "level2 = " << level_func1 << ";" << std::endl;
      output << this->indent[3] << "}" << std::endl;
      // Replace alias names
      level_func1 = std::string("level");
      level_func2 = std::string("level2");
      index_func1 = std::string("index");
      index_func2 = std::string("index2");
    }
    // Loop to evaluate the base function on the left, right or mid points of the other basefunction
    for (int i = 0; i < 1 + static_cast<int>(use_implicit_zero); i++) {
      if (use_level_cache) {
        // Reuse h values from host
        output << this->indent[3] << "h = hs[" << level_func2 << "];" << std::endl;
      } else {
        // Calculate h
        output << this->indent[3] << "h = 1.0 / (1 << " << level_func2 << ");" << std::endl;
      }
      // Calculate u
      output << this->indent[3] << "u = (1 << " << level_func1 << ");" << std::endl;
      // Check whether we will just need to calculate umid, or umid uright and uleft
      if (use_less) {
        // Calculate just umid
        output << this->indent[3] << "umid = u * h * (" << index_func2 << ") - " << index_func1
               << ";" << std::endl;
        output << this->indent[3] << "umid = fabs(umid);" << std::endl;
        output << this->indent[3] << "umid = 1.0-umid;" << std::endl;
        if (!use_fabs_instead_of_fmax)
          output << this->indent[3] << "umid = fmax(umid,0.0);" << std::endl;
        else
          output << this->indent[3] << "umid = (umid + fabs(umid));" << std::endl;
        // Add integral to result sum
        if (i == 0)
          output << this->indent[3] << "sum = h*(umid);" << std::endl;
        else
          output << this->indent[3] << "sum += h*(umid);" << std::endl;
      } else {
        // Calculate umid, uright and uleft
        output << this->indent[3] << "umid = u * h * (" << index_func2 << ") - " << index_func1
               << ";" << std::endl;
        output << this->indent[3] << "umid = fabs(umid);" << std::endl;
        output << this->indent[3] << "umid = 1.0-umid;" << std::endl;
        if (!use_fabs_instead_of_fmax)
          output << this->indent[3] << "umid = fmax(umid,0.0);" << std::endl;
        else
          output << this->indent[3] << "umid = (umid + fabs(umid));" << std::endl;
        output << this->indent[3] << "uright = u*h*(" << index_func2 << " + 1) - " << index_func1
               << ";" << std::endl;
        output << this->indent[3] << "uright = fabs(uright);" << std::endl;
        output << this->indent[3] << "uright = 1.0-uright;" << std::endl;
        if (!use_fabs_instead_of_fmax)
          output << this->indent[3] << "uright = fmax(uright,0.0);" << std::endl;
        else
          output << this->indent[3] << "uright = (uright + fabs(uright));" << std::endl;
        output << this->indent[3] << "uleft = u*h*(" << index_func2 << " - 1) - " << index_func1
               << ";" << std::endl;
        output << this->indent[3] << "uleft = fabs(uleft);" << std::endl;
        output << this->indent[3] << "uleft = 1.0-uleft;" << std::endl;
        if (!use_fabs_instead_of_fmax)
          output << this->indent[3] << "uleft = fmax(uleft,0.0);" << std::endl;
        else
          output << this->indent[3] << "uleft = (uleft + fabs(uleft));" << std::endl;
        // Add integral to result sum
        if (i == 0)
          output << this->indent[3] << "sum = h/3.0*(umid + uleft + uright);" << std::endl;
        else
          output << this->indent[3] << "sum += h/3.0*(umid + uleft + uright);" << std::endl;
      }
      // Swap aliases
      level_func2 = std::string("point_level_block") + std::to_string(block) + std::string("[dim]");
      level_func1 = std::string("starting_points[i* ") + std::to_string(dimensions) +
                    std::string("*2+2*dim+1]");
      index_func2 =
          std::string("point_indices_block") + std::to_string(block) + std::string("[dim]");
      index_func1 = std::string("starting_points[i* ") + std::to_string(dimensions) +
                    std::string("*2+2*dim]");
      if (useLocalMemory) {
        level_func1 =
            std::string("level_local[i* ") + std::to_string(dimensions) + std::string("+dim]");
        index_func1 =
            std::string("indices_local[i* ") + std::to_string(dimensions) + std::string("+dim]");
      }
    }
    // Check whether we need to do something about base functions with the same level
    if (use_less) {
      if (!do_not_use_ternary) {
        if (use_implicit_zero) {
          // Use ternary operator to multiply with 1/3
          output << this->indent[3] << "sum *= " << level_func2 << " == " << level_func1
                 << " ? 1.0/3.0 : 1.0;" << std::endl;
        } else {
          // Use ternary operator to multiply with 2/3
          output << this->indent[3] << "sum *= " << level_func2 << " == " << level_func1
                 << " ? 2.0/3.0 : 1.0;" << std::endl;
        }
      } else {
        // decrement counter of same levels by one if the two levels do not match
        output << this->indent[3] << "same_levels -= min((int)(abs(level_local[i*" << dimensions
               << "+dim] - point_level_block" << block << "[dim])),(int)(1));" << std::endl;
      }
    } else if (!use_implicit_zero) {
      if (!do_not_use_ternary) {
        output << this->indent[3] << "sum *= " << level_func2 << " == " << level_func1
               << " ? 2.0 : 1.0;" << std::endl;
      } else {
        // decrement counter of same levels by one if the two levels do not match
        output << this->indent[3] << "same_levels -= min((int)(abs(level_local[i*" << dimensions
               << "+dim] - point_level_block" << block << "[dim])),(int)(1));" << std::endl;
      }
    }
    // Update cell integral
    output << this->indent[3] << "zellenintegral*=sum;" << std::endl;
    output << this->indent[2] << "}" << std::endl;
    // Update cell integral with missing factors
    if (do_not_use_ternary) {
      if (use_less) {
        output << this->indent[2] << "zellenintegral *= divisors[same_levels];" << std::endl;
        if (!use_implicit_zero) {
          output << this->indent[2] << "zellenintegral *= (1 << same_levels);" << std::endl;
        }
      } else if (!use_implicit_zero) {
        output << this->indent[2] << "zellenintegral *= (1 << same_levels);" << std::endl;
      }
    }
    // Mutliply with corresponding alpha value
    if (useLocalMemory) {
      output << this->indent[2] << "gesamtint_block" << block
             << " += zellenintegral*alpha_local[i];" << std::endl;
    } else {
      output << this->indent[2] << "gesamtint_block" << block << " += zellenintegral*alpha[i];"
             << std::endl;
    }
    return output.str();
  }

 public:
  explicit SourceBuilderMult(json::Node &kernelConfiguration)
      : kernelConfiguration(kernelConfiguration),
        dataBlockSize(1),
        use_level_cache(false),
        use_less(true),
        do_not_use_ternary(false),
        use_implicit_zero(true),
        use_fabs_instead_of_fmax(false),
        preprocess_positions(false),
        unroll_dim(false) {
    if (kernelConfiguration.contains("LOCAL_SIZE"))
      localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    if (kernelConfiguration.contains("KERNEL_USE_LOCAL_MEMORY"))
      useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    if (kernelConfiguration.contains("KERNEL_DATA_BLOCKING_SIZE"))
      dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    if (kernelConfiguration.contains("USE_LEVEL_CACHE"))
      use_level_cache = kernelConfiguration["USE_LEVEL_CACHE"].getBool();
    if (kernelConfiguration.contains("USE_LESS_OPERATIONS"))
      use_less = kernelConfiguration["USE_LESS_OPERATIONS"].getBool();
    if (kernelConfiguration.contains("DO_NOT_USE_TERNARY"))
      do_not_use_ternary = kernelConfiguration["DO_NOT_USE_TERNARY"].getBool();
    if (kernelConfiguration.contains("USE_IMPLICIT"))
      use_implicit_zero = kernelConfiguration["USE_IMPLICIT"].getBool();
    if (kernelConfiguration.contains("USE_FABS"))
      use_fabs_instead_of_fmax = kernelConfiguration["USE_FABS"].getBool();
    if (kernelConfiguration.contains("PREPROCESS_POSITIONS")) {
      preprocess_positions = kernelConfiguration["PREPROCESS_POSITIONS"].getBool();
      // These two options are not compatible
      if (preprocess_positions) use_level_cache = false;
    }
    if (kernelConfiguration.contains("UNROLL_DIM")) {
      unroll_dim = kernelConfiguration["UNROLL_DIM"].getBool();
    }
  }

  /// Generates the opencl source code for the density matrix-vector multiplication
  std::string generateSource(size_t dimensions, size_t problemsize) {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("DensityMultiplication.cl");
    }

    std::stringstream sourceStream;

    if (this->floatType().compare("double") == 0) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    if (!preprocess_positions) {
      sourceStream << "void multdensity(__global const int *starting_points,";
    } else {
      sourceStream << "void multdensity(__global const int *hs_inverses, __global const "
                   << this->floatType() << " *hs," << std::endl
                   << " __global const " << this->floatType() << " *positions,";
    }
    sourceStream << "__global const " << this->floatType() << " *alpha, __global "
                 << this->floatType() << " *result,const " << this->floatType()
                 << " lambda, const int startid";
    if (use_level_cache) sourceStream << ", __global " << this->floatType() << " *hs";
    if (do_not_use_ternary) sourceStream << ", __global " << this->floatType() << " *divisors";
    sourceStream << ")" << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << this->indent[0] << "int gridindex = startid + get_global_id(0);" << std::endl
                 << this->indent[0] << "__private int local_id = get_local_id(0);" << std::endl;
    sourceStream << save_from_global_to_private(dimensions);
    sourceStream << this->indent[0] << "__private int teiler = 0;" << std::endl;
    sourceStream << this->indent[0] << "__private " << this->floatType() << " h = 1.0 / 3.0;"
                 << std::endl;
    sourceStream << this->indent[0] << "__private " << this->floatType() << " umid = 0.0;"
                 << std::endl;
    if (!use_less) {
      sourceStream << this->indent[0] << "__private " << this->floatType() << " uright = 0.0;"
                   << std::endl;
      sourceStream << this->indent[0] << "__private " << this->floatType() << " uleft = 0.0;"
                   << std::endl;
    }
    sourceStream << this->indent[0] << "__private " << this->floatType() << " sum = 0.0;"
                 << std::endl;
    sourceStream << this->indent[0] << "__private int u= 0;" << std::endl;
    for (size_t block = 0; block < dataBlockSize; block++)
      sourceStream << this->indent[0] << this->floatType() << " gesamtint_block" << block
                   << " = 0.0;" << std::endl;
    // Store points in local memory
    if (useLocalMemory && !preprocess_positions) {
      sourceStream << this->indent[0] << "__local "
                   << "int indices_local[" << localWorkgroupSize * dimensions << "];" << std::endl
                   << this->indent[0] << "__local "
                   << "int level_local[" << localWorkgroupSize * dimensions << "];" << std::endl
                   << this->indent[0] << "__local " << this->floatType() << " alpha_local["
                   << localWorkgroupSize << "];" << std::endl
                   << this->indent[0] << "for (int group = 0; group < "
                   << problemsize / localWorkgroupSize << "; group++) {" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int j = 0; j <     " << dimensions << " ; j++) {"
                   << std::endl
                   << this->indent[2] << "indices_local[local_id * " << dimensions
                   << " + j] = starting_points[group * " << localWorkgroupSize * dimensions * 2
                   << "  + local_id * " << dimensions * 2 << " + 2 * j];" << std::endl
                   << this->indent[2] << "level_local[local_id * " << dimensions
                   << " + j] = starting_points[group * " << localWorkgroupSize * dimensions * 2
                   << "  + local_id * " << dimensions * 2 << " + 2 * j + 1];" << std::endl
                   << this->indent[1] << "}" << std::endl
                   << this->indent[1] << "alpha_local[local_id] = alpha[group * "
                   << localWorkgroupSize << "  + local_id ];" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int i = 0 ; i < " << localWorkgroupSize << "; i++) {"
                   << std::endl
                   << this->indent[2] << "__private " << this->floatType();
    } else if (preprocess_positions && !unroll_dim) {
      // declare local arrays for grid point positions, and hs and hs inverses
      sourceStream << this->indent[0] << "__local " << this->floatType() << " positions_local["
                   << localWorkgroupSize * dimensions << "];" << std::endl
                   << this->indent[0] << "__local " << this->floatType() << " hs_local["
                   << localWorkgroupSize * dimensions << "];" << std::endl
                   << this->indent[0] << "__local int hinverses_local["
                   << localWorkgroupSize * dimensions << "];" << std::endl
                   << this->indent[0] << "__local " << this->floatType() << " alpha_local["
                   << localWorkgroupSize << "];" << std::endl;
      // start loop
      sourceStream << this->indent[0] << "for (int group = 0; group < "
                   << problemsize / localWorkgroupSize << "; group++) {" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int j = 0; j <     " << dimensions << " ; j++) {"
                   << std::endl;
      // get hs inverse
      sourceStream << this->indent[2] << "hinverses_local[local_id * " << dimensions
                   << " + j] = hs_inverses[group * " << localWorkgroupSize * dimensions
                   << " + local_id*" << dimensions << " + j];" << std::endl;
      // get hs
      sourceStream << this->indent[2] << "hs_local[local_id * " << dimensions
                   << " + j] = hs[group * " << localWorkgroupSize * dimensions << " + local_id*"
                   << dimensions << " + j];" << std::endl;
      // get positions
      sourceStream << this->indent[2] << "positions_local[local_id * " << dimensions
                   << " + j] = positions[group * " << localWorkgroupSize * dimensions
                   << " + local_id*" << dimensions << " + j];" << std::endl;
      sourceStream << "}" << this->indent[1] << "alpha_local[local_id] = alpha[group * "
                   << localWorkgroupSize << "  + local_id ];" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int i = 0 ; i < " << localWorkgroupSize << "; i++) {"
                   << std::endl
                   << this->indent[2] << "__private " << this->floatType();
    } else if (preprocess_positions && unroll_dim) {
      for (size_t dim = 0; dim < dimensions; ++dim) {
        sourceStream << this->indent[0] << "__local " << this->floatType() << " positions_local_dim"
                     << dim << "[" << localWorkgroupSize << "];" << std::endl
                     << this->indent[0] << "__local " << this->floatType() << " hs_local"
                     << "_dim" << dim << "[" << localWorkgroupSize << "];" << std::endl
                     << this->indent[0] << "__local int hinverses_local"
                     << "_dim" << dim << "[" << localWorkgroupSize << "];" << std::endl;
      }
      sourceStream << this->indent[0] << "__local " << this->floatType() << " alpha_local["
                   << localWorkgroupSize << "];" << std::endl;
      // start loop
      sourceStream << this->indent[0] << "for (int group = 0; group < "
                   << problemsize / localWorkgroupSize << "; group++) {" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      for (size_t dim = 0; dim < dimensions; ++dim) {
        // get hs inverse
        sourceStream << this->indent[2] << "hinverses_local_dim" << dim
                     << "[local_id] = hs_inverses[group * " << localWorkgroupSize * dimensions
                     << " + local_id*" << dimensions << " + " << dim << "];" << std::endl;
        // get hs
        sourceStream << this->indent[2] << "hs_local_dim" << dim << "[local_id] = hs[group * "
                     << localWorkgroupSize * dimensions << " + local_id*" << dimensions << " + "
                     << dim << "];" << std::endl;
        // get positions
        sourceStream << this->indent[2] << "positions_local_dim" << dim
                     << "[local_id] = positions[group * " << localWorkgroupSize * dimensions
                     << " + local_id*" << dimensions << " + " << dim << "];" << std::endl;
      }
      sourceStream << this->indent[1] << "alpha_local[local_id] = alpha[group * "
                   << localWorkgroupSize << "  + local_id ];" << std::endl
                   << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << this->indent[1] << "for (int i = 0 ; i < " << localWorkgroupSize << "; i++) {"
                   << std::endl
                   << this->indent[2] << "__private " << this->floatType();
    } else {
      sourceStream << this->indent[1] << "for (int i = 0 ; i < " << problemsize << "; i++) {"
                   << std::endl
                   << this->indent[2] << "__private " << this->floatType();
    }

    // Generate body for each element in the block
    for (size_t block = 0; block < dataBlockSize; block++) {
      if (preprocess_positions) {
        sourceStream << " zellenintegral = 1.0;" << std::endl;
        if (!unroll_dim) {
          // Loop over all dimensions
          sourceStream << this->indent[2] << "for(private int dim = 0;dim< " << dimensions
                       << ";dim++) {" << std::endl;
          // Calculate distance between grid points in this dimensions
          sourceStream << this->indent[2] << "__private " << this->floatType()
                       << " distance = fabs(point_positions_block" << block << "[dim] - "
                       << "positions_local[i * " << dimensions << " + dim]);" << std::endl;
          // Calculate first integral (level_point < level_local)
          sourceStream << this->indent[2] << "sum = 1.0 - distance * point_hinverses_block" << block
                       << "[dim]; " << std::endl;
          sourceStream << this->indent[2] << "sum *= hs_local[i *" << dimensions << " + dim]; "
                       << std::endl;
          sourceStream << this->indent[2] << "sum = max(sum, (" << this->floatType() << ") 0.0); "
                       << std::endl;
          // Calculate second integral (level_point > level_local)
          sourceStream << this->indent[2] << "sum += max((" << this->floatType()
                       << ")(point_hs_block" << block << "[dim] * (1.0 - hinverses_local[i * "
                       << dimensions << " + dim] * distance)), (" << this->floatType() << ")0.0);"
                       << std::endl;
          // Update cell integral
          if (this->floatType().compare("double") == 0) {
            sourceStream << this->indent[3] << "zellenintegral*=sum*select((" << this->floatType()
                         << " )1.0,"
                         << " (" << this->floatType() << " )h, (long)(point_hinverses_block"
                         << block << "[dim] == hinverses_local[i * " << dimensions << " + dim]));"
                         << std::endl;
          } else {
            sourceStream << this->indent[3] << "zellenintegral*=sum*select((" << this->floatType()
                         << " )1.0,"
                         << " (" << this->floatType() << " )h, (uint)(point_hinverses_block"
                         << block << "[dim] == hinverses_local[i * " << dimensions << " + dim]));"
                         << std::endl;
          }
          sourceStream << this->indent[2] << "}" << std::endl;
        } else {
          for (size_t dim = 0; dim < dimensions; ++dim) {
            // Calculate distance between grid points in this dimensions
            sourceStream << this->indent[2] << "__private " << this->floatType() << " distance_dim"
                         << dim << " = fabs(point_positions_block" << block << "[" << dim << "] - "
                         << "positions_local_dim" << dim << "[i]);" << std::endl;
            // Calculate first integral (level_point < level_local)
            sourceStream << this->indent[2] << "sum = 1.0 - distance_dim" << dim
                         << " * point_hinverses_block" << block << "[" << dim << "]; " << std::endl;
            sourceStream << this->indent[2] << "sum *= hs_local_dim" << dim << "[i]; " << std::endl;
            sourceStream << this->indent[2] << "sum = max(sum, (" << this->floatType() << ")0.0); "
                         << std::endl;
            // Calculate second integral (level_point > level_local)
            sourceStream << this->indent[2] << "sum += max((" << this->floatType()
                         << ")(point_hs_block" << block << "[" << dim
                         << "] * (1.0 - hinverses_local_dim" << dim << "[i] * distance_dim" << dim
                         << ")), (" << this->floatType() << ")0.0);" << std::endl;
            // Update cell integral
            if (this->floatType().compare("double") == 0) {
              sourceStream << this->indent[3] << "zellenintegral*=sum*select((" << this->floatType()
                           << " )1.0,"
                           << " (" << this->floatType() << " )h, (long)(point_hinverses_block"
                           << "[" << dim << "] == hinverses_local_dim" << dim << "[i]));"
                           << std::endl;
            } else {
              sourceStream << this->indent[3] << "zellenintegral*=sum*select((" << this->floatType()
                           << " )1.0,"
                           << " (" << this->floatType() << " )h, (uint)(point_hinverses_block"
                           << "[" << dim << "] == hinverses_local_dim" << dim << "[i]));"
                           << std::endl;
            }
          }
        }
        sourceStream << this->indent[2] << "gesamtint_block" << block
                     << " += zellenintegral*alpha_local[i];" << std::endl;
      } else {
        sourceStream << calculate_matrix_entry(block, dimensions) << std::endl;
      }
    }
    // Close group loop
    if (useLocalMemory) sourceStream << this->indent[1] << "}" << std::endl;
    sourceStream << this->indent[0] << "}" << std::endl;
    for (size_t block = 0; block < dataBlockSize; ++block) {
      if (!use_fabs_instead_of_fmax || preprocess_positions)
        sourceStream << this->indent[0] << "result[get_global_id(0) * " << dataBlockSize << " + "
                     << block << "] = gesamtint_block" << block << ";" << std::endl;
      else
        sourceStream << this->indent[0] << "result[get_global_id(0) * " << dataBlockSize << " + "
                     << block << "] = gesamtint_block" << block << " / " << (1 << dimensions) << ";"
                     << std::endl;
      sourceStream << this->indent[0] << "result[get_global_id(0) * " << dataBlockSize << " + "
                   << block << "] += alpha[gridindex * " << dataBlockSize << " + " << block << "]*"
                   << "lambda;" << std::endl;
    }

    sourceStream << "}" << std::endl;

    if (kernelConfiguration.contains("WRITE_SOURCE")) {
      if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
        this->writeSource("DensityMultiplication.cl", sourceStream.str());
      }
    }
    return sourceStream.str();
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
