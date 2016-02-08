/*
 * hierarchization.hpp
 *
 *  Created on: 02.07.2014
 *      Author: P. Butz
 */

#ifndef HIERARCHIZATION_HPP_
#define HIERARCHIZATION_HPP_

#include "boost/lexical_cast.hpp"
#include <cstdlib>
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/combigrid_ultils.hpp"
#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"

/*
 * Instead of having private static functions, I put these functions in an
 * unnamed namespace. So, they are not accessible from outside the file, as well.
 * In the general case, this would have the advantage, that we can change
 * the declaration of these functions without changing the declaration of the
 * class. Thus, we avoid recompilation of all files that use the class.
 * However, as here everything is in the header file, this probably won't
 * have any effect. At least for building with make, I don't know for scons...
 */
namespace {

using namespace combigrid;

template<typename FG_ELEMENT>
inline static void hierarchize1DUnoptimizedNoBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim);

template<typename FG_ELEMENT>
inline static void hierarchize1DUnoptimizedBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim);

template<typename FG_ELEMENT>
inline void dehierarchize1DUnoptimizedBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim);

template<typename FG_ELEMENT>
inline void dehierarchize1DUnoptimizedNoBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim);

template<typename FG_ELEMENT>
static void printValues2DArr(FG_ELEMENT* val, IndexType size, IndexType offset,
    IndexType n0);

template<typename FG_ELEMENT>
static void printValues(combigrid::FullGrid<FG_ELEMENT>& fg);
}

namespace combigrid {

class Hierarchization {

public:
  // inplace hierarchization
  template<typename FG_ELEMENT>
  static void hierarchize(FullGrid<FG_ELEMENT>& fg) {
    assert(!fg.isHierarchized());

    const DimType d = fg.getDimension();
    // number of grid points in each dimension
    IndexVector n(d);
    for (DimType i = 0; i < fg.getDimension(); ++i) {
      if (fg.returnBoundaryFlags()[i] == false) {
        n[i] = static_cast<IndexType>(std::pow(2, fg.getLevels()[i])) - 1;
      } else {
        n[i] = static_cast<IndexType>(std::pow(2, fg.getLevels()[i])) + 1;
      }
    }

    const IndexType size = fg.getNrElements(); // number of grid points

    IndexType stride = 1;
    IndexType jump;
    lldiv_t divresult;

    theStatsContainer()->setTimerStart("hierarchize_dim_0");

    //   dimension 1 separate as start of each pole is easier to calculate
    IndexType ndim = n[0];
    IndexType nbrOfPoles = size / ndim;
    IndexType start = -100;
    if (fg.returnBoundaryFlags()[0] == false) {
#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
      for (IndexType kk = 0; kk < nbrOfPoles; kk++) {
        start = kk * ndim;
        hierarchize1DUnoptimizedNoBoundary(fg, start, 1, ndim, 0);
      }
    } else {
#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
      for (IndexType kk = 0; kk < nbrOfPoles; kk++) {
        start = kk * ndim;
        hierarchize1DUnoptimizedBoundary(fg, start, 1, ndim, 0);
      }
    }
    // end dimension 1
    theStatsContainer()->setTimerStop("hierarchize_dim_0");

    for (DimType dim = 1; dim < d; dim++) { // hierarchize for all dims
      theStatsContainer()->setTimerStart(
          "hierarchize_dim_" + boost::lexical_cast<std::string>(dim));

      stride *= ndim;
      ndim = n[dim];
      jump = stride * ndim;
      nbrOfPoles = size / ndim;
//      std::cout << "dim " << dim << "size " << size << "ndim " << ndim << "nbrpoles" << nbrOfPoles << std::endl;
      if (fg.returnBoundaryFlags()[dim] == false) {
#pragma omp parallel for schedule(static) firstprivate (divresult, stride, dim, ndim, nbrOfPoles, start, jump)
        for (IndexType nn = 0; nn < nbrOfPoles; nn++) { // integer operations form bottleneck here -- nested loops are twice as slow
          divresult = std::lldiv(nn, stride);
          start = divresult.quot * jump + divresult.rem;
          hierarchize1DUnoptimizedNoBoundary(fg, start, stride, ndim, dim);
        }
      } else {
#pragma omp parallel for schedule(static) firstprivate (divresult, stride, dim, ndim, nbrOfPoles, start, jump)
        for (IndexType nn = 0; nn < nbrOfPoles; nn++) { // integer operations form bottleneck here -- nested loops are twice as slow
          divresult = std::lldiv(nn, stride);
          start = divresult.quot * jump + divresult.rem;
          hierarchize1DUnoptimizedBoundary(fg, start, stride, ndim, dim);
        }
      }

      theStatsContainer()->setTimerStop(
          "hierarchize_dim_" + boost::lexical_cast<std::string>(dim));

    } // end loop over dimension 2 to d

    //printValues(fg);

    // set fg to hierarchized
    fg.isHierarchized_ = true;

    return;
  }

  // inplace dehierarchization
  template<typename FG_ELEMENT>
  static void dehierarchize(FullGrid<FG_ELEMENT>& fg) {
    assert(fg.isHierarchized());

    const DimType d = fg.getDimension();
    // number of grid points in each dimension
    IndexVector n(d);
    for (DimType i = 0; i < fg.getDimension(); ++i) {
      if (fg.returnBoundaryFlags()[i] == false) {
        n[i] = static_cast<IndexType>(std::pow(2, fg.getLevels()[i])) - 1;
      } else {
        n[i] = static_cast<IndexType>(std::pow(2, fg.getLevels()[i])) + 1;
      }
    }

    const IndexType size = fg.getNrElements(); // number of grid points

    IndexType stride = 1;
    IndexType jump;
    lldiv_t divresult;

    //   dimension 1 separate as start of each pole is easier to calculate
    IndexType ndim = n[0];
    IndexType nbrOfPoles = size / ndim;
    IndexType start = -100;
    if (fg.returnBoundaryFlags()[0] == false) {
#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
      for (IndexType kk = 0; kk < nbrOfPoles; kk++) {
        start = kk * ndim;
        dehierarchize1DUnoptimizedNoBoundary(fg, start, 1, ndim, 0);
      }
    } else {
#pragma omp parallel for schedule(static) firstprivate (ndim, nbrOfPoles, start)
      for (IndexType kk = 0; kk < nbrOfPoles; kk++) {
        start = kk * ndim;
        dehierarchize1DUnoptimizedBoundary(fg, start, 1, ndim, 0);
      }
    }
    // end dimension 1
//	    printValues(fg);

    for (DimType dim = 1; dim < d; dim++) { // hierarchize for all dims
      stride *= ndim;
      ndim = n[dim];
      jump = stride * ndim;
      nbrOfPoles = size / ndim;
      if (fg.returnBoundaryFlags()[dim] == false) {
#pragma omp parallel for schedule(static) firstprivate (divresult, stride, dim, ndim, nbrOfPoles, start, jump)
        for (IndexType nn = 0; nn < nbrOfPoles; nn++) { // integer operations form bottleneck here -- nested loops are twice as slow
          divresult = std::lldiv(nn, stride);
          start = divresult.quot * jump + divresult.rem;
          dehierarchize1DUnoptimizedNoBoundary(fg, start, stride, ndim, dim);
        }
      } else {
#pragma omp parallel for schedule(static) firstprivate (divresult, stride, dim, ndim, nbrOfPoles, start, jump)
        for (IndexType nn = 0; nn < nbrOfPoles; nn++) { // integer operations form bottleneck here -- nested loops are twice as slow
          divresult = std::lldiv(nn, stride);
          start = divresult.quot * jump + divresult.rem;
          dehierarchize1DUnoptimizedBoundary(fg, start, stride, ndim, dim);
        }
      }
//	      printValues(fg);

    } // end loop over dimension 2 to d

    // set fg to unhierarchized
    fg.isHierarchized_ = false;

    return;

  }

};

} // namespace combigrid

namespace {

using namespace combigrid;
/**
 * Solves the Hierarchization for a specific 1D problem
 */
template<typename FG_ELEMENT>
inline void hierarchize1DUnoptimizedNoBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim) {

  IndexType steps;
  IndexType ctr;
  IndexType offset, parentOffset;
  IndexType stepsize;
  IndexType parOffsetStrided;

  FG_ELEMENT* val = fg.getData();

  // ssa variables
  FG_ELEMENT val1 = 0, val2 = 0, val3 = 0, parL = 0, parR = 0;

  IndexType ll = fg.getLevels()[dim];
  steps = (1 << (ll - 1));
  offset = 0;
  stepsize = 2;
  parentOffset = 1;

  for (ll--; ll > 0; ll--) {
    parOffsetStrided = parentOffset * stride;
    val[start + offset * stride] -= 0.5
        * val[start + offset * stride + parOffsetStrided];
    offset += stepsize;
    parL = 0.5 * val[start + offset * stride - parOffsetStrided];
    for (ctr = 1; ctr < steps - 1; ctr++) {
      val1 = val[start + offset * stride];
      parR = 0.5 * val[start + offset * stride + parOffsetStrided];
      val2 = val1 - parL;
      val3 = val2 - parR;
      val[start + offset * stride] = val3;
      parL = parR;
      offset += stepsize;
    }
    val[start + offset * stride] -= parL;
    steps = steps >> 1;
    offset = (1 << (fg.getLevels()[dim] - ll)) - 1;
    parentOffset = stepsize;
    stepsize = stepsize << 1;
  }
  return;
}

template<typename FG_ELEMENT>
inline void hierarchize1DUnoptimizedBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim) {

  IndexType steps;
  IndexType ctr;
  IndexType offset, parentOffset;
  IndexType stepsize;
  IndexType parOffsetStrided;

  FG_ELEMENT* val = fg.getData();

  // ssa variables
  FG_ELEMENT val1 = -100, val2 = -100, val3 = -100, parL = -100, parR = -100;

  IndexType ll = fg.getLevels()[dim];
  steps = (1 << (ll - 1));
  offset = 1; // 1 da boundary
  stepsize = 2;
  parentOffset = 1;

  for (ll--; ll > -1; ll--) { // hier is index um 1 geschiftet da vorher level 2 manuell behandelt wurde
    parOffsetStrided = parentOffset * stride;
    parL = 0.5 * val[start + offset * stride - parOffsetStrided];
    for (ctr = 0; ctr < steps; ctr++) {
      val1 = val[start + offset * stride];
      parR = 0.5 * val[start + offset * stride + parOffsetStrided];
      val2 = val1 - parL;
      val3 = val2 - parR;
      val[start + offset * stride] = val3;
      parL = parR;
      offset += stepsize;
    }
    steps = steps >> 1;
    offset = (1 << (fg.getLevels()[dim] - ll)); // boundary case
    parentOffset = stepsize;
    stepsize = stepsize << 1;
  }
  return;
}

template<typename FG_ELEMENT>
inline void dehierarchize1DUnoptimizedNoBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim) {

//	      std::cout << "start " << start
//	                << "stride " << stride
//	                << "size " << size
//	                << "dim " << dim << std::endl;

  IndexType ll;
  IndexType steps;
  IndexType ctr;
  IndexType offset, parentOffset;
  IndexType stepsize;
  IndexType parOffsetStrided;

  FG_ELEMENT* val = fg.getData();

  // ssa variables
  FG_ELEMENT val1 = -100, val2 = -100, val3 = -100, parL = -100, parR = -100;

  LevelType maxL = fg.getLevels()[dim];
//				int maxL = l[dim];
//			start = start +stride; // nun mit offset = 1;start war vorher der erste gitterpunkt (welcher randpunkt ist) und ist nun der erste zu hierarchisierende Gitterpunkt.
  steps = 2;
  offset = (1 << (maxL - 2)) - 1; // offset =1 da boundary.
  stepsize = (1 << (maxL - 1));
  parentOffset = (1 << (maxL - 2));

  for (IndexType ll = 2; ll <= maxL; ll++) {
    parOffsetStrided = parentOffset * stride;
    parL = 0.5 * val[start + offset * stride + parOffsetStrided]; //eigentlich rechts hier aber einmal umkopieren sparen.
    val[start + offset * stride] += parL;
    offset += stepsize;
//				parL= 0.5*val[start+offset*stride -parOffsetStrided];

    for (ctr = 1; ctr < steps - 1; ctr++) {
      val1 = val[start + offset * stride];
      parR = 0.5 * val[start + offset * stride + parOffsetStrided];
      val2 = val1 + parL;
      val3 = val2 + parR;
      val[start + offset * stride] = val3;
      parL = parR;
      offset += stepsize;
    }
    val[start + offset * stride] += parL;
    steps = steps << 1;
    offset = (1 << (maxL - (ll + 1))) - 1; // boundary case
    parentOffset = parentOffset >> 1;
    stepsize = stepsize >> 1;
//				std::cout << "after hierarchizing level " << ll << std::endl;
//				printValues(fg);
  }
  return;
}

template<typename FG_ELEMENT>
inline void dehierarchize1DUnoptimizedBoundary(
    combigrid::FullGrid<FG_ELEMENT>& fg, IndexType start, IndexType stride,
    IndexType size, DimType dim) {
  IndexType ll;
  IndexType steps;
  IndexType ctr;
  IndexType offset, parentOffset;
  IndexType stepsize;
  IndexType parOffsetStrided;

  FG_ELEMENT* val = fg.getData();

  // ssa variables
  FG_ELEMENT val1 = -100, val2 = -100, val3 = -100, parL = -100, parR = -100;

  LevelType maxL = fg.getLevels()[dim];
//				int maxL = l[dim];
//			start = start +stride; // nun mit offset = 1;start war vorher der erste gitterpunkt (welcher randpunkt ist) und ist nun der erste zu hierarchisierende Gitterpunkt.
  steps = 1;
  offset = (1 << (maxL - 1)); // offset =1 da boundary.
  stepsize = (1 << maxL);
  parentOffset = (1 << (maxL - 1));

  for (LevelType ll = 1; ll <= maxL; ll++) { // couting with offset of one as level 2 was hierarchized manually before
    // just convers setting of strides and co.
    parOffsetStrided = parentOffset * stride;
    parL = 0.5 * val[start + offset * stride - parOffsetStrided];
    ;
    for (ctr = 0; ctr < steps; ctr++) {
      val1 = val[start + offset * stride];
      parR = 0.5 * val[start + offset * stride + parOffsetStrided];
      val2 = val1 + parL;
      val3 = val2 + parR;
      val[start + offset * stride] = val3;
      parL = parR;
      offset += stepsize;
    }
//				val[start+offset*stride] -= parR;
    steps = steps << 1;
    offset = (1 << (maxL - (ll + 1))); // boundary case
    parentOffset = parentOffset >> 1;
    stepsize = stepsize >> 1;
  }
  return;
}

template<typename FG_ELEMENT>
static void printValues2DArr(FG_ELEMENT* val, IndexType size, IndexType offset,
    IndexType n0) {
  for (IndexType i = 1; i <= size; i++) {
    std::cout << val[offset + i - 1] << "\t";
    if (i % n0 == 0)
      std::cout << std::endl;
  }
  std::cout << std::endl;
  return;
}

template<typename FG_ELEMENT>
static void printValues(combigrid::FullGrid<FG_ELEMENT>& fg) {
  printf("\n");
  const DimType d = fg.getDimension();
  // number of grid points in each dimension
  const IndexVector& n = fg.getSizes();

  FG_ELEMENT* val = fg.getData();

  const IndexType size = fg.getNrElements();

  if (d <= 2)
    printValues2DArr(val, size, 0, n[0]);
  else {
    IndexType * currentLevels = new IndexType[d];
    IndexType chunkSize = n[0] * n[1];
    IndexType nbrOfChunks = size / chunkSize;
    for (IndexType ctr = 0; ctr < nbrOfChunks; ctr++) {
      std::cout << std::endl;
      printValues2DArr(val, chunkSize, ctr * chunkSize, n[0]);
      // rest is only for formatting - make an additioanl new line when increasing the level for d >= 2
      currentLevels[2]++;
      for (IndexType dd = 2; dd < d; dd++) {
        if (currentLevels[dd] == n[dd]) {
          currentLevels[dd] = 0;
          currentLevels[dd + 1]++;
          std::cout << std::endl;
        }
      }
    }
  }
  return;
}

} /* namespace */

#endif /* HIERARCHISERIELL_HPP_ */
