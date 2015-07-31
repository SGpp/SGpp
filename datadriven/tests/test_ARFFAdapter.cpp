/*#define BOOST_TEST_DYN_LINK

#include <unistd.h>
#include <limits.h>

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

using namespace SGPP::base;

std::string get_selfpath() {
  char buff[PATH_MAX];
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff) - 1);

  if (len != -1) {
    buff[len] = '\0';
    return std::string(buff);
  } else {
    return "/fail/";
  }
}

BOOST_AUTO_TEST_SUITE(testARFFAdapter)

BOOST_AUTO_TEST_CASE(testLoadData) {
SGPP::float_t testPoints[10][3] = {
  {0.307143, 0.130137, 0.050000},
  {0.365584, 0.105479, 0.050000},
  {0.178571, 0.201027, 0.050000},
  {0.272078, 0.145548, 0.050000},
  {0.318831, 0.065411, 0.050000},
  {0.190260, 0.086986, 0.050000},
  {0.190260, 0.062329, 0.072500},
  {0.120130, 0.068493, 0.072500},
  {0.225325, 0.056164, 0.072500},
  {0.213636, 0.050000, 0.072500}
};

SGPP::float_t testValues[10] = { -1., 1., 1., 1., 1., 1., -1., -1., -1., -1.};

  //std::string filename = get_selfpath() + "/datasets/liver-disorders_normalized.arff.gz";
  SGPP::datadriven::Dataset readData = SGPP::datadriven::ARFFTools::readARFF("datadriven/tests/datasets/liver-disorders_normalized.arff.gz");
  DataVector* classes = readData.getClasses();
  DataMatrix* train = readData.getTrainingData();
  size_t size = classes->getSize();
  size_t dim = train->getNrows();

  //DataVector testVector = DataVector(dim);
  for (size_t rowIdx = 0; rowIdx < size; rowIdx++) {
    //train->getRow(rowIdx, testVector);
    for (size_t colIdx = 0; colIdx < dim; colIdx++) {
      BOOST_CHECK_EQUAL(train->get(rowIdx, colIdx), testPoints[rowIdx][colIdx]);
    }

    BOOST_CHECK_EQUAL(classes->get(rowIdx), testValues[rowIdx]);
  }
}


BOOST_AUTO_TEST_SUITE_END()
*/