#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include "test_datadrivenCommon.hpp"

using namespace SGPP::base;
using namespace SGPP::datadriven;


DataMatrix* generateBTMatrix(Grid* grid, DataMatrix& training) {

  GridStorage* storage = grid->getStorage();

  OperationMultipleEval* b = SGPP::op_factory::createOperationMultipleEval(*grid,
                             training);

  DataVector alpha(storage->size());
  DataVector temp(training.getNrows());

  // create BT matrix
  DataMatrix* m = new DataMatrix(training.getNrows(), storage->size());

  for (size_t i = 0; i < storage->size(); i++) {
    temp.setAll(0.0);
    alpha.setAll(0.0);
    alpha[i] = 1.0;
    b->mult(alpha, temp);

    m->setColumn(i, temp);

  }

  return m;
}

void compareBTMatrices(DataMatrix* m1, DataMatrix* m2) {

#if USE_DOUBLE_PRECISION
  double tolerance = 1E-5;
#else
  double tolerance = 1E-2;
#endif

  // check dimensions
  BOOST_CHECK_EQUAL(m1->getNrows(), m2->getNrows());
  BOOST_CHECK_EQUAL(m1->getNcols(), m2->getNcols());

  size_t rows = m1->getNrows(); //was n

  size_t cols = m1->getNcols(); //was m

  // check row sum
  DataVector v(cols);

  std::vector<SGPP::float_t> values;

  for (size_t i = 0; i < rows; i++) {
    m1->getRow(i, v);
    values.push_back(v.sum());

    //    std::cout << "row[" << i << "]: ";
    //    for (size_t j = 0; j < v.getSize(); j++) {
    //      if (j > 0) {
    //        std::cout << ", ";
    //      }
    //      std::cout << "(" << v[j] << ")";
    //    }
    //    std::cout << std::endl;
  }

  //  std::cout << "------------------------------" << std::endl;

  //std::sort(values.begin(), values.end());

  std::vector<SGPP::float_t> valuesReference;

  for (size_t i = 0; i < rows; i++) {
    m2->getRow(i, v);
    valuesReference.push_back(v.sum());

    //    std::cout << "row[" << i << "]: ";
    //    for (size_t j = 0; j < v.getSize(); j++) {
    //      if (j > 0) {
    //        std::cout << ", ";
    //      }
    //      std::cout << "(" << v[j] << ")";
    //    }
    //    std::cout << std::endl;
  }

  //std::sort(valuesReference.begin(), valuesReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_CLOSE(values[i], valuesReference[i], tolerance);
  }

  // check col sum
  values.clear();

  DataVector vRows(rows);

  for (size_t i = 0; i < cols; i++) {
    m1->getColumn(i, vRows);
    values.push_back(vRows.sum());


  }

  std::sort(values.begin(), values.end());

  valuesReference.clear();

  for (size_t i = 0; i < cols; i++) {
    m2->getColumn(i, vRows);
    valuesReference.push_back(vRows.sum());


  }

  std::sort(valuesReference.begin(), valuesReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_CLOSE(values[i], valuesReference[i], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE(TestOperationBTModLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_ausgeklappt_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {

  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {

  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_ausgeklappt_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}


BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {

  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {

  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinearBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}


BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {

  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {

  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinearTruncatedBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {

  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}


BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {

  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {

  size_t level = 2;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {

  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string
  referenceMatrixFileName("datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix &trainingData = dataset.getTrainingData();

  size_t dim = dataset.getDimension();

  Grid* grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()
