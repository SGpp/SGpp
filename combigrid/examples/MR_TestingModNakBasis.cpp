//#include <sgpp/base/datatypes/DataVector.hpp>
//#include <sgpp/base/grid/Grid.hpp>
//#include <sgpp/combigrid/GeneralFunction.hpp>
//#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
//#include <sgpp/combigrid/pce/HierarchicalStochasticCollocation.hpp>
//#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
//#include <sgpp/optimization/sle/solver/Auto.hpp>
//#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
//#include <sgpp/optimization/tools/Printer.hpp>
//
//#include <iostream>
//#include <memory>
//#include "../../base/src/sgpp/base/grid/type/NakBsplineModifiedGrid.hpp"
//
// double l2ErrorWithBoundary(std::shared_ptr<sgpp::base::Grid> surrogateGrid,
//                           sgpp::base::DataVector alpha, double
//                           (*objFunc)(sgpp::base::DataVector), size_t dim) {
//  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*surrogateGrid, alpha);
//  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
//  sgpp::base::GridStorage& gridStorage = grid->getStorage();
//  size_t level = 9;
//  grid->getGenerator().full(level);
//  sgpp::base::DataVector diffSquare(gridStorage.getSize(), 0.0);
//  for (size_t i = 0; i < gridStorage.getSize(); i++) {
//    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
//    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
//    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
//      p[j] = gp.getStandardCoordinate(j);
//    }
//    diffSquare[i] = std::pow(fabs(objFunc(p) - sparseGridSurrogate.eval(p)), 2);
//  }
//  return sqrt(diffSquare.sum());
//}
//
// double l2ErrorNoBoundary(std::shared_ptr<sgpp::base::Grid> surrogateGrid,
//                         sgpp::base::DataVector alpha, double (*objFunc)(sgpp::base::DataVector),
//                         size_t dim) {
//  sgpp::optimization::InterpolantScalarFunction sparseGridSurrogate(*surrogateGrid, alpha);
//  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
//  sgpp::base::GridStorage& gridStorage = grid->getStorage();
//  size_t level = 9;
//  grid->getGenerator().full(level);
//  sgpp::base::DataVector diffSquare(gridStorage.getSize(), 0.0);
//  for (size_t i = 0; i < gridStorage.getSize(); i++) {
//    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
//    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
//    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
//      p[j] = gp.getStandardCoordinate(j);
//    }
//    diffSquare[i] = std::pow(fabs(objFunc(p) - sparseGridSurrogate.eval(p)), 2);
//  }
//  return sqrt(diffSquare.sum());
//}
//
//// double objectiveFunction(sgpp::base::DataVector v) { return sin(v[0] * v[1]) * cos(v[0] *
///v[1]); / }
// double objectiveFunction(sgpp::base::DataVector v) { return 1; }
// double weightFunction(double x) { return 1; }
//
// int main() {
//  size_t dim = 1;
//  size_t degree = 3;
//  sgpp::combigrid::WeightFunctionsCollection weights(
//      dim, sgpp::combigrid::SingleFunction(weightFunction));
//  sgpp::base::DataVector bounds(2 * dim, 0.0);
//  for (size_t d = 0; d < dim; d++) {
//    bounds[2 * d] = 0.0;
//    bounds[2 * d + 1] = 1.0;
//  }
//  for (size_t level = 1; level < 6; level++) {
//    std::shared_ptr<sgpp::base::Grid> ModNakGrid(
//        sgpp::base::Grid::createNakBsplineModifiedGrid(dim, degree));
//    ModNakGrid->getGenerator().regular(level);
//    sgpp::base::GridStorage& gridStorage = ModNakGrid->getStorage();
//    sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
//    for (size_t i = 0; i < gridStorage.getSize(); i++) {
//      sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
//      sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
//      for (size_t j = 0; j < gridStorage.getDimension(); j++) {
//        p[j] = gp.getStandardCoordinate(j);
//      }
//      f_values[i] = objectiveFunction(p);
//    }
//
//    sgpp::optimization::sle_solver::Auto sleSolver;
//    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
//    sgpp::optimization::HierarchisationSLE hierSLE(*ModNakGrid);
//    sgpp::base::DataVector alpha(ModNakGrid->getSize());
//    if (!sleSolver.solve(hierSLE, f_values, alpha)) {
//      std::cout << "Solving failed!" << std::endl;
//    }
//    sgpp::optimization::InterpolantScalarFunction modGridSurrogate(*ModNakGrid, alpha);
//
//    //----------------
//    std::shared_ptr<sgpp::base::Grid> BoundaryNakGrid(
//        sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, degree));
//    BoundaryNakGrid->getGenerator().regular(level);
//    sgpp::base::GridStorage& boundaryGridStorage = BoundaryNakGrid->getStorage();
//    sgpp::base::DataVector boundary_f_values(boundaryGridStorage.getSize(), 0.0);
//    for (size_t i = 0; i < boundaryGridStorage.getSize(); i++) {
//      sgpp::base::GridPoint& gp = boundaryGridStorage.getPoint(i);
//      sgpp::base::DataVector p(boundaryGridStorage.getDimension(), 0.0);
//      for (size_t j = 0; j < boundaryGridStorage.getDimension(); j++) {
//        p[j] = gp.getStandardCoordinate(j);
//      }
//      boundary_f_values[i] = objectiveFunction(p);
//    }
//
//    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
//    sgpp::optimization::HierarchisationSLE boundarySLE(*BoundaryNakGrid);
//    sgpp::base::DataVector boundaryAlpha(BoundaryNakGrid->getSize());
//    if (!sleSolver.solve(boundarySLE, boundary_f_values, boundaryAlpha)) {
//      std::cout << "Solving failed!" << std::endl;
//    }
//    sgpp::optimization::InterpolantScalarFunction boundaryGridSurrogate(*BoundaryNakGrid,
//                                                                        boundaryAlpha);
//    //--------------
//
//    //    double modL2errorB = l2ErrorWithBoundary(ModNakGrid, alpha, *objectiveFunction, dim);
//    //    double boundaryL2errorB =
//    //        l2ErrorWithBoundary(BoundaryNakGrid, boundaryAlpha, *objectiveFunction, dim);
//    //    double modL2error = l2ErrorNoBoundary(ModNakGrid, alpha, *objectiveFunction, dim);
//    //    double boundaryL2error =
//    //        l2ErrorNoBoundary(BoundaryNakGrid, boundaryAlpha, *objectiveFunction, dim);
//
//    //    std::cout << "level: " << level << std::endl;
//    //    std::cout << "GP | mod: " << ModNakGrid->getSize()
//    //              << " boundary: " << BoundaryNakGrid->getSize() << std::endl;
//    //    std::cout << "Error B| mod: " << modL2errorB << " boundary: " << boundaryL2errorB <<
//    //    std::endl; std::cout << "Error| mod: " << modL2error << " boundary: " << boundaryL2error
//    //    << std::endl;
//
//    //    sgpp::combigrid::HierarchicalStochasticCollocation hsc(BoundaryNakGrid, boundaryAlpha,
//    //    weights,
//    //                                                           bounds);
//  }
//  //-----------------------------------------------
//  sgpp::combigrid::HierarchicalStochasticCollocation hsc(
//      sgpp::base::GridType::NakBsplineModified, dim,
//      sgpp::combigrid::MultiFunction(objectiveFunction), weights, bounds, degree);
//  hsc.refineRegular(1);
//
//  std::cout << "coeff: " << hsc.getCoefficients().toString() << std::endl;
//
//  std::cout << "#gp: " << hsc.numGridPoints() << std::endl;
//  std::cout << "mean: " << hsc.mean() << std::endl;
//  std::cout << "var: " << hsc.variance() << std::endl;
//  std::cout << "\n";
//
//  hsc.refineSurplusAdaptive(3);
//
//  std::cout << "coeff: " << hsc.getCoefficients().toString() << std::endl;
//
//  std::cout << "#gp: " << hsc.numGridPoints() << std::endl;
//  std::cout << "mean: " << hsc.mean() << std::endl;
//  std::cout << "var: " << hsc.variance() << std::endl;
//  std::cout << "\n";
//  return 0;
//}
