// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/GridPrinterForStretching.hpp>
#include <sgpp/base/exception/tool_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>


namespace sgpp {
namespace base {

GridPrinterForStretching::GridPrinterForStretching(Grid& SparseGrid) :
  GridPrinter(SparseGrid) {
  myGrid = &SparseGrid;
}

GridPrinterForStretching::~GridPrinterForStretching() {
}

void GridPrinterForStretching::printGridDomain(DataVector& alpha,
    std::string tFilename, BoundingBox& GridArea, size_t PointsPerDimension) {
  throw tool_exception("GridPrinterForStretching::printGridDomain : "
                           "This printer does not support BoundingBox, "
                           "use printGridDomainStretching instead!");
}

void GridPrinterForStretching::printGridDomainStretching(DataVector& alpha,
    std::string tFilename, Stretching& GridArea, size_t PointsPerDimension) {
  DimensionBoundary dimOne;
  DimensionBoundary dimTwo;
  std::ofstream fileout;

  if (myGrid->getSize() > 0) {
    if (myGrid->getDimension() != 2) {
      throw tool_exception(
        "GridPrinterForStretching::printGridDomainStretching : "
        "The grid has more not two dimensions. Thus it cannot be printed!");
    } else {
      // Open filehandle
      fileout.open(tFilename.c_str());
      std::unique_ptr<OperationEval> myEval = sgpp::op_factory::createOperationEval(*myGrid);

      dimOne = GridArea.getBoundary(0);
      dimTwo = GridArea.getBoundary(1);

      for (double i = dimOne.leftBoundary; i <= dimOne.rightBoundary;
           i += ((dimOne.rightBoundary - dimOne.leftBoundary) /
                 static_cast<double>
                 (PointsPerDimension))) {
        for (double j = dimTwo.leftBoundary; j <= dimTwo.rightBoundary;
             j += ((dimTwo.rightBoundary - dimTwo.leftBoundary) /
                   static_cast<double>
                   (PointsPerDimension))) {
          std::vector<double> point;
          point.push_back(i);
          point.push_back(j);
          fileout << i << " " << j << " " << myEval->eval(alpha, point) <<
                  std::endl;
        }

        fileout << std::endl;
      }

      // close filehandle
      fileout.close();
    }
  } else {
    throw tool_exception(
      "GridPrinterForStretching::printGridDomainStretching : "
      "The grid has no dimensions. Thus it cannot be printed!");
  }
}

void GridPrinterForStretching::printGrid(DataVector& alpha,
    std::string tFilename, size_t PointsPerDimension) {
  DimensionBoundary dimOne;
  DimensionBoundary dimTwo;
  std::ofstream fileout;

  if (myGrid->getSize() > 0) {
    if (myGrid->getDimension() > 2) {
      throw tool_exception("GridPrinter::printGrid : "
                               "The grid has more than two dimensions. "
                               "Thus it cannot be printed!");
    } else {
      // Open filehandle
      fileout.open(tFilename.c_str());
      std::unique_ptr<OperationEval> myEval = sgpp::op_factory::createOperationEval(*myGrid);

      if (myGrid->getDimension() == 1) {
        dimOne = myGrid->getStretching().getBoundary(0);

        double inc_x = ((dimOne.rightBoundary - dimOne.leftBoundary) /
                         (static_cast<double>(PointsPerDimension) - 1.0));

        size_t points = PointsPerDimension;

        for (size_t i = 0; i < points; i++) {
          std::vector<double> point;
          point.push_back(((static_cast<double>(i))*inc_x));
          fileout << ((static_cast<double>(i))*inc_x) << " " << myEval->eval(alpha,
                  point) << std::endl;
        }
      } else if (myGrid->getDimension() == 2) {
        dimOne = myGrid->getStretching().getBoundary(0);
        dimTwo = myGrid->getStretching().getBoundary(1);

        double offset_x = dimOne.leftBoundary;
        double offset_y = dimTwo.leftBoundary;
        double inc_x = ((dimOne.rightBoundary - dimOne.leftBoundary) /
                         (static_cast<double>(PointsPerDimension) - 1.0));
        double inc_y = ((dimTwo.rightBoundary - dimTwo.leftBoundary) /
                         (static_cast<double>(PointsPerDimension) - 1.0));

        size_t points = PointsPerDimension;

        for (size_t i = 0; i < points; i++) {
          for (size_t j = 0; j < points; j++) {
            std::vector<double> point;
            point.push_back(offset_x + ((static_cast<double>(i))*inc_x));
            point.push_back(offset_y + ((static_cast<double>(j))*inc_y));
            fileout << (offset_x + (static_cast<double>(i))*inc_x) << " " <<
                    (offset_y + (static_cast<double>(j))*inc_y) << " " <<
                    myEval->eval(alpha, point) << std::endl;
          }

          fileout << std::endl;
        }
      }

      // close filehandle
      fileout.close();
    }

  } else {
    throw tool_exception("GridPrinterForStretching::printGrid : "
                             "The grid has no dimensions. "
                             "Thus it cannot be printed!");
  }
}

void GridPrinterForStretching::printSparseGrid(DataVector& alpha,
    std::string tFilename, bool bSurplus) {
  DataVector temp(alpha);
  double tmp = 0.0;
  size_t dim = myGrid->getDimension();
  std::ofstream fileout;

  // Do Dehierarchisation, is specified
  if (bSurplus == false) {
    sgpp::op_factory::createOperationHierarchisation(*myGrid)->doDehierarchisation(temp);
  }

  // Open filehandle
  fileout.open(tFilename.c_str());

  for (size_t i = 0; i < myGrid->getSize(); i++) {
    std::string coords =
      myGrid->getStorage().get(i)->getCoordsStringStretching(
        myGrid->getStretching());
    std::stringstream coordsStream(coords);

    for (size_t j = 0; j < dim; j++) {
      coordsStream >> tmp;
      fileout << tmp << " ";
    }

    fileout << temp[i] << std::endl;
  }

  fileout.close();
}

void GridPrinterForStretching::printSparseGridExpTransform(DataVector& alpha,
    std::string tFilename, bool bSurplus) {
  DataVector temp(alpha);
  double tmp = 0.0;
  size_t dim = myGrid->getDimension();
  std::ofstream fileout;

  // Do Dehierarchisation, is specified
  if (bSurplus == false) {
    sgpp::op_factory::createOperationHierarchisation(*myGrid)->doDehierarchisation(temp);
  }

  // Open filehandle
  fileout.open(tFilename.c_str());

  for (size_t i = 0; i < myGrid->getSize(); i++) {
    std::string coords =
      myGrid->getStorage().get(i)->getCoordsStringStretching(
        myGrid->getStretching());
    std::stringstream coordsStream(coords);

    for (size_t j = 0; j < dim; j++) {
      coordsStream >> tmp;
      fileout << exp(tmp) << " ";
    }

    fileout << temp[i] << std::endl;
  }

  fileout.close();
}

}  // namespace base
}  // namespace sgpp
