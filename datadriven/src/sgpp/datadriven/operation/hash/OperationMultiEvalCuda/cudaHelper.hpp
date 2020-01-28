// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

#ifndef CUDA_ERROR_CHECK
#define CUDA_ERROR_CHECK

///@cond DOXY_IGNORE // NOLINT()
namespace sgpp {
namespace datadriven {
namespace OpMultiEvalCudaDetail {

/// Wrapper for __cudaSafeCall
#define CudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
/// Wrapper for __cudaCheckError
#define CudaCheckError() __cudaCheckError(__FILE__, __LINE__)

/// Error catcher for CUDA API calls
inline void __cudaSafeCall(cudaError err, const char* file, const int line) {
#ifdef CUDA_ERROR_CHECK
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }
#endif
  return;
}

/// Checker for previously occured errors in GPU kernel code
inline void __cudaCheckError(const char* file, const int line) {
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }
  // More careful checking. However, this will affect performance.
  // Comment away if needed.
  err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n", file, line,
            cudaGetErrorString(err));
    exit(-1);
  }
#endif
  return;
}

/// Template class managing host and device memory
template <typename _dtype>
class HostDevPtr {
 public:
  /// Constructor does not allocate any memory
  HostDevPtr() {
    _host = nullptr;
    _dev = nullptr;
    _size = 0;
    _ref = false;
  }

  /// Destructor frees memory on host and device if allocated
  ~HostDevPtr() {
    if (_host && !_ref) delete[](_host);
    if (_dev) cudaFree(_dev);
  }

  /// Frees memory if allocated and sets size to 0
  void Free() {
    if (_host && !_ref) {
      delete[](_host);
      _host = nullptr;
      _size = 0;
      _ref = false;
    }
    if (_dev) {
      cudaFree(_dev);
      _dev = nullptr;
    }
  }

  /// Copies data from host to device
  void CopyToDev() {
    if (_host && _dev) cudaMemcpy(_dev, _host, sizeof(_dtype) * _size, cudaMemcpyHostToDevice);
  }
  /// Copies data from device to host
  void CopyToHost() {
    if (_host && _dev) cudaMemcpy(_host, _dev, sizeof(_dtype) * _size, cudaMemcpyDeviceToHost);
  }
  /// Allocates an array on the host
  void HostAlloc(const size_t& size) {
    if (size == _size) return;
    if (_host && !_ref) delete[](_host);
    if (_dev) {
      cudaFree(_dev);
      _dev = nullptr;
    }
    _host = new _dtype[size];
    _size = size;
    _ref = false;
  }
  /// Sets extern array as source. Does not allocate new memory
  void HostRef(_dtype* ptr, const size_t& size) {
    if (_host) delete[](_host);
    if (_dev) {
      cudaFree(_dev);
      _dev = nullptr;
    }
    _host = ptr;
    _size = size;
    _ref = true;
  }
  /// Allocates the corresponding amount of memory on the device
  void DevAlloc() {
    if (!_host || _dev) return;
    cudaMalloc(reinterpret_cast<void**>(&_dev), sizeof(_dtype) * _size);
    CudaCheckError();
  }
  /// Excesses elements on the host
  _dtype& operator[](const size_t& idx) {
    if (!_host || idx >= _size) throw;
    return _host[idx];
  }
  /// Returns the number of elements
  const size_t& size() const { return _size; }
  /// Returns the pointer to the host memory
  _dtype* host() const { return _host; }
  /// Returns the pointer to the device memory
  _dtype* dev() const { return _dev; }
  /// Writes the host data to a binary file
  void Save(const char* file) {
    FILE* handle = fopen(file, "wb");
    fwrite(_host, _size, sizeof(_dtype), handle);
    fclose(handle);
  }
  /// Loads data from a binary file
  void Load(const char* file) {
    FILE* handle = fopen(file, "rb");
    fpos_t p;
    fseek(handle, 0, SEEK_END);
    fgetpos(handle, &p);
    Free();
    _size = p.__pos / sizeof(_dtype);
    _ref = false;
    _host = new _dtype[_size];
    rewind(handle);
    if (fread(_host, _size, sizeof(_dtype), handle) == _size) printf("OK\n");
    fclose(handle);
  }

 private:
  _dtype* _host;
  _dtype* _dev;
  size_t _size;
  bool _ref;
};

}  // namespace OpMultiEvalCudaDetail
}  // namespace datadriven
}  // namespace sgpp
///@endcond // NOLINT()

#endif  // CUDA_ERROR_CHECK
