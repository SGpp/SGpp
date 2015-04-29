// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OCLKERNELBASE_HPP
#define OCLKERNELBASE_HPP

#include <CL/cl.h>
#include <string>
#include <sgpp/parallel/datadriven/basis/common/KernelBase.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    template<typename T> struct getType {};

    template<> struct getType<float> {
      static std::string asString() {
        return "float";
      }
      static std::string constSuffix() {
        return "f";
      }
      static std::string intAsString() {
        return "uint";
      }
    };

    template<> struct getType<double> {
      static std::string asString() {
        return "double";
      }
      static std::string constSuffix() {
        return "";
      }
      static std::string intAsString() {
        return "ulong";
      }
    };

    class OCLKernelBase {
      private:
        virtual std::string generateSourceMultTrans(size_t dims, size_t local_workgroup_size) = 0;
        virtual std::string generateSourceMult(size_t dims, size_t local_workgroup_size) = 0;
      public:
        cl_int createMultTrans(size_t dims, size_t local_workgroup_size, cl_context context, size_t num_devices, cl_device_id* device_ids, cl_kernel* kernel) {
          std::string program_src = generateSourceMultTrans(dims, local_workgroup_size);
          return OCLKernelBase::buildKernel(program_src, "multTransOCL", context, num_devices, device_ids, kernel);
        }

        cl_int createMult(size_t dims, size_t local_workgroup_size, cl_context context, size_t num_devices, cl_device_id* device_ids, cl_kernel* kernel) {
          std::string program_src = generateSourceMult(dims, local_workgroup_size);
          return OCLKernelBase::buildKernel(program_src, "multOCL", context, num_devices, device_ids, kernel);
        }

      private:

        /**
         * @brief buildKernel builds the program that is represented by @a program_src and creates @a num_devices kernel objects
         * that are stored into the array @a kernel (must be already allocated with at least @a num_devices )
         *
         * @param program_src the source of the program to compile
         * @param kernel_name name of the kernel function (in program_src) to create the kernel for
         * @param context OpenCL context
         * @param num_devices number of OpenCL devices
         * @param device_ids array with device ids, necessary for displaying build info
         * @param kernel already allocated array: the resulting kernels are put into this array, one for each device (=> at least num_devices entries)
         * @return
         */
        static inline cl_int buildKernel(const std::string& program_src, const char* kernel_name,
                                         cl_context context, size_t num_devices,
                                         cl_device_id* device_ids, cl_kernel* kernel) {
          cl_int err;

          //std::cout << program_src << std::endl;

          // setting the program
          const char* kernel_src = program_src.c_str();
          cl_program program = clCreateProgramWithSource(context, 1, &kernel_src, NULL, &err);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: Failed to create program! Error Code: " << err << std::endl;
            throw SGPP::base::operation_exception("OCL Error: Failed to create program!");
          }

          std::string build_opts;
#ifndef NO_OCL_OPTS
          build_opts = "-cl-finite-math-only -cl-fast-relaxed-math ";
#else
          build_opts = "-cl-opt-disable ";
#endif

          // compiling the program
          err = clBuildProgram(program, 0, NULL, build_opts.c_str(), NULL, NULL);

          if (err != CL_SUCCESS) {
            std::cout << "OCL Error: OpenCL Build Error. Error Code: " << err << std::endl;

            size_t len;
            char buffer[4096];

            // get the build log
            clGetProgramBuildInfo(program, device_ids[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
            std::cout << "--- Build Log ---" << std::endl << buffer << std::endl;
            return err;
          }

          // creating the kernel
          for (size_t i = 0; i < num_devices; i++) {
            kernel[i] = clCreateKernel(program, kernel_name, &err);

            if (err != CL_SUCCESS) {
              std::cout << "OCL Error: Failed to create kernel! Error Code: " << err << std::endl;
              return err;
            }
          }

          if (program) {
            clReleaseProgram(program);
          }

          return CL_SUCCESS;
        }
    };
  }
}
#endif // OCLKERNELBASE_HPP
