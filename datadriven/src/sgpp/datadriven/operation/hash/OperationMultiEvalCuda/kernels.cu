#include <stdint.h>
#include "kernels.cuh"

/// Device variables for magic numbers
__constant__ double __idxtable[4];
__constant__ int32_t __idxtablei[4];

__constant__ uint32_t __dDIM[1];

#define __DIM __dDIM[0]

/// Shared gridnode objects. Works up to Dim = 32. Be aware of limited shared memory on GPU before increasing values!
__shared__ gridnode_t snode[32];

/// Initializes magic numbers and dimension constants on the GPU
void initCudaConstants (const uint32_t DIM) {
  cudaMemcpyToSymbol(__idxtablei, idxtablei, sizeof(int32_t)*4);
  cudaMemcpyToSymbol(__idxtable, idxtable, sizeof(double)*4);
  cudaMemcpyToSymbol(__dDIM, &DIM, sizeof(uint32_t));
}

/// Evaluation of a DIM-dim polynomial basis function
/// The node of the basis function is passed by the global shared snode
/// For further code comments see evalBasis_1
__device__ double evalBasis (double*& pos, uint32_t start) {
  double eval = 1.0;
  double x;
  uint32_t id;
  int32_t root;
  for (uint32_t d=0;d<__DIM;++d) {
    x = pos[start + d] * snode[d].level2;
    id = x;
    id |= 1;
    if (id != snode[d].index) return 0.0;
    x = double(id + 1) - x;
    eval *= x;
    root = -1;
    x -= 2.0;
    for (int32_t j=2;j<snode[d].grad;j<<=1) {
      eval *= x/root;
      root += (__idxtablei[id&3] * j);
			x    += (__idxtablei[id&3] * j);
			id>>=1;
    }
    eval *= x/root;
  }
  return eval;
}

/// Evaluation of (DIM-1)-dim polynomial basis function. It skipps the first dimension
/// For further code comments see evalBasis_1
__device__ double evalBasis_d (double* pos, gridnode_t* node) {
  double eval = 1.0;
  double x;
  uint32_t id;
  int32_t root;
  for (uint32_t d=1;d<__DIM;++d) {
    x = pos[d] * node[d].level2;
    id = x;
    id |= 1;
    x = double(id + 1) - x;
    eval *= x;
    root = -1;
    x -= 2.0;
    for (int32_t j=2;j<node[d].grad;j<<=1) {
      eval *= x/root;
      root += (__idxtablei[id&3] * j);
			x    += (__idxtablei[id&3] * j);
			id>>=1;
    }
    eval *= x/root;
  }
  return eval;
}

/// Evaluation of a one-dim polynomial basis function.
__device__ double evalBasis_1 (double pos, gridnode_t& node) {
  // This function iterativley constructs the Legandre polynomial
  // The localy nearest ancestors of the node are used as roots
  double eval = 1.0;
  uint32_t id;
  int32_t root;
  // Scale position to level
  pos *= node.level2;
  // Compute index according to the eval position
  id = pos;
  id |= 1;
  // Go to right support boundary as 1st root
  pos = double(id + 1) - pos;
  eval *= pos;
  root = -1;
  // Go to left support boundary for 2nd root
  pos -= 2.0;
  for (int32_t j=2;j<node.grad;j<<=1) {
    // Multiply new root
    eval *= pos/root;
    // Use magic numbers and shifted level to find next ancestor
    root += (__idxtablei[id&3] * j);
    pos    += (__idxtablei[id&3] * j);
    // Shift the level
    id>>=1;
  }
  // Multiply last root
  eval *= pos/root;
  return eval;
}

__global__ void gpu_zindex (double* pos, uint64_t* index) {
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
  pos = &pos[idx*__DIM];
  uint64_t val = 0;
  //uint32_t bits = 64/__DIM;
	for (uint32_t i=0;i<64/__DIM;++i) {
		for (uint32_t d=0;d<__DIM;++d) {
			val <<= 1;
			if (pos[d] >= 1.0) val |= 1;
			else
				val |= uint64_t(pos[d]*(1<<i))&1;
		}
	}
	index[idx] = val;
	//index[idx] = zvalue(&pos[idx*__DIM]);
}

/// Standard evaluation. Each parallel thread computes one evaluation.
__global__ void gpu_eval (double *res, gridnode_t* node, double* a, double* pos,
    uint32_t* limit, uint32_t subcnt, uint32_t* subs) {
  uint32_t idx = threadIdx.x + blockIdx.x*blockDim.x;
  res = &res[idx];
  pos = &pos[idx*__DIM];
  limit = &limit[idx*__DIM];
  
  double tmp, poly;
  double data = 0.0;
  uint32_t path;
  
  // Iterate through all possible subspaces with levelsum = n + d - 1
  // The subspaces are pre-computed
  for (uint32_t j=0;j<subcnt*__DIM;j+=__DIM) {
    // Start with root node
    idx = GRID_START;
    
    // Follow path in all but thre first dimension
    for (uint32_t d=1;d<__DIM;++d) {
      // Skip if evaluation point is outside of any support
      if (subs[j + d] > limit[d]) {
        idx = GRID_END;
        break;
      }
      // Follow path to children
      path = pos[d] * (1<<subs[j + d]);
      path = __brev(path) >> (32 - subs[j + d]);
      for (uint32_t l=1;l<subs[j + d];++l) {
        if (idx == GRID_END) break;
        idx = node[idx*__DIM + d].child[path&1];
        path >>= 1;
      }
      // Skip if leaf is reached
      if (idx == GRID_END) break;
    }
    if (idx == GRID_END) continue;
    
    // Do the actual evaluation along the remaining one-dim subtree
    path = pos[0] * (1<<subs[j]);
    path = __brev(path) >> (32 - subs[j]);
    tmp = 0.0;
    poly = evalBasis_d(pos, &node[idx*__DIM]);
    // Go down tree and accumulate weighted basis functions
    for (uint32_t l=1;l<=subs[j];++l) {
			if (idx == GRID_END) break;
			tmp += a[idx] * evalBasis_1(pos[0],node[idx*__DIM]);
			idx = node[idx*__DIM].child[(path&1)];
			path >>= 1;
		}
		// Geather for all possible subspaces
		data += tmp * poly;
  }
  res[0] = data;
}

/// Transposed evaluation with optimized streaming apporach and additional FMA. Each thread block computes one result entry
__global__ void gpu_transevel (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit, double* b, double c, double M_) {
  __shared__ double val[CUDA_BLOCKSIZE];
  __shared__ limit_t slimit;
  double tmp = 0.0;
  if (threadIdx.x < __DIM) {
    snode[threadIdx.x] = node[blockIdx.x*__DIM + threadIdx.x];
  }
  slimit = limit[blockIdx.x];
  __syncthreads();
  
  // Streaming loop strided within the thread block
  // The upper and lower limits can be pre-computed if the dataset is aligned along a space-filling curve
  // Otherwise use 0 and M as limits
  for (uint32_t i=slimit.lower + threadIdx.x; i < slimit.upper;
      i += CUDA_BLOCKSIZE) {
    tmp += y[i] * evalBasis(pos, i*__DIM);
  }
  val[threadIdx.x] = tmp;
  
  // Fan-in all sub results
  __syncthreads();
	#if CUDA_BLOCKSIZE > 512
    if (threadIdx.x < 512)
		val[threadIdx.x] += val[threadIdx.x + 512];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 256
	if (threadIdx.x < 256)
		val[threadIdx.x] += val[threadIdx.x + 256];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 128
	if (threadIdx.x < 128)
		val[threadIdx.x] += val[threadIdx.x + 128];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 64
	if (threadIdx.x < 64)
		val[threadIdx.x] += val[threadIdx.x + 64];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 32
	if (threadIdx.x < 32)
		val[threadIdx.x] += val[threadIdx.x + 32];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 16];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 8];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 4];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 2];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 1];
	__syncthreads();
	
	// Compute final result
	if (threadIdx.x == 0) {
    tmp = val[0] * M_;
		a[blockIdx.x] = (c*b[blockIdx.x]) + tmp;
  }
}

/// Compute subspace limits for each evaluation point
__global__ void gpu_preprocess (uint32_t* limit, gridnode_t* node, double* pos,
    uint32_t maxlevel) {
  uint32_t idx = threadIdx.x + blockIdx.x*blockDim.x;
  pos = &pos[idx*__DIM];
  limit = &limit[idx*__DIM];
  uint32_t path;
  for (uint32_t d=0;d<__DIM;++d) {
    idx = GRID_START;
    path = pos[d]*(1<<maxlevel);
		path = __brev(path) >> (32 - maxlevel);
    limit[d] = maxlevel;
    for (uint32_t l=1;l<=maxlevel;++l) {
      if (idx == GRID_END) {
          limit[d] = l;
          break;
      }
      idx = node[idx*__DIM + d].child[path&1];
      path >>= 1;
    }
  }
}

/// Transposed evaluation with optimized streaming approach
/// See the other gpu_transeval kernel for code comments
__global__ void gpu_transevel (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit) {
  __shared__ double val[CUDA_BLOCKSIZE];
  //__shared__ gridnode_t snode[32];
  double tmp = 0.0;
  if (threadIdx.x < __DIM) {
    snode[threadIdx.x] = node[blockIdx.x*__DIM + threadIdx.x];
  }
  __syncthreads();
  for (uint32_t i=limit[blockIdx.x].lower + threadIdx.x; i < limit[blockIdx.x].upper;
      i += CUDA_BLOCKSIZE) {
    //tmp += y[i] * evalBasis(&pos[i*__DIM], snode);//&node[blockIdx.x*__DIM]);
    tmp += y[i] * evalBasis(pos, i*__DIM);
  }
  val[threadIdx.x] = tmp;
  __syncthreads();
	#if CUDA_BLOCKSIZE > 512
    if (threadIdx.x < 512)
		val[threadIdx.x] += val[threadIdx.x + 512];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 256
	if (threadIdx.x < 256)
		val[threadIdx.x] += val[threadIdx.x + 256];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 128
	if (threadIdx.x < 128)
		val[threadIdx.x] += val[threadIdx.x + 128];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 64
	if (threadIdx.x < 64)
		val[threadIdx.x] += val[threadIdx.x + 64];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	#if CUDA_BLOCKSIZE > 32
	if (threadIdx.x < 32)
		val[threadIdx.x] += val[threadIdx.x + 32];
	__syncthreads();
	#endif // CUDA_BLOCKSIZE
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 16];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 8];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 4];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 2];
	__syncthreads();
	if (threadIdx.x < 16)
		val[threadIdx.x] += val[threadIdx.x + 1];
	__syncthreads();
	if (threadIdx.x == 0)
		a[blockIdx.x] = val[threadIdx.x];
}

/// Compute streaming boundary limits for datasets aligned along a Morton order curve
__global__ void gpu_zbound (uint64_t* idx_p, gridnode_t* node, limit_t* limit, uint32_t M) {
  // Use the orderd dataset as binary search tree to find position of upper and lower
  // support corner within the space-filling curve. Each thread computes the limits for
  // one grid node.
  uint32_t idx = threadIdx.x + blockIdx.x*blockDim.x;
  node = &node[idx*__DIM];
  uint64_t idx_l, idx_h;
  double xl, xh;
  idx_h = 0;
  idx_l = 0;
  uint32_t bits = 64/__DIM;
	for (uint32_t i=0;i<bits;++i) {
		for (uint32_t d=0;d<__DIM;++d) {
      xl = node[d].x - 1.0/node[d].level2;
      xh = node[d].x + 0.8/node[d].level2;
			idx_h <<= 1;
      idx_h |= uint64_t(xh*(1<<i))&1;
      idx_l <<= 1;
      idx_l |= uint64_t(xl*(1<<i))&1;
		}
	}
  uint32_t l,h;
	l = 0;
	h = 0;
  for (uint32_t i=M/2;i>0;i>>=1) {
		if (idx_p[i+l] <  idx_l)
			l += i;
		if (idx_p[i+h] <= idx_h)
				h += i;
	}
	limit[idx].lower  = l&(0xFFFFFFFF-(CUDA_BLOCKSIZE-1));
	limit[idx].upper = min((h+2)|(CUDA_BLOCKSIZE-1),M);
}
