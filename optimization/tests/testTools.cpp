#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cstdio>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/optimization/tools/Math.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include "ObjectiveFunctions.hpp"

const bool use_double_precision =
#if USE_DOUBLE_PRECISION
  true;
#else
  false;
#endif /* USE_DOUBLE_PRECISION */

using namespace SGPP;
using namespace SGPP::optimization;

void gridEqualityTest(base::Grid& grid1, base::Grid& grid2) {
  base::GridStorage& storage1 = *grid1.getStorage();
  base::GridStorage& storage2 = *grid2.getStorage();
  const size_t d = storage1.dim();
  BOOST_CHECK_EQUAL(d, storage2.dim());
  const size_t n = storage1.size();
  BOOST_CHECK_EQUAL(n, storage2.size());

  for (size_t k = 0; k < n; k++) {
    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(storage1[k]->getLevel(t),
                        storage2[k]->getLevel(t));
      BOOST_CHECK_EQUAL(storage1[k]->getIndex(t),
                        storage2[k]->getIndex(t));
    }
  }
}

SGPP::float_t calculateMean(std::vector<SGPP::float_t>& x) {
  SGPP::float_t mean = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    mean += x[i];
  }

  mean /= static_cast<SGPP::float_t>(x.size());
  return mean;
}

SGPP::float_t calculateVariance(std::vector<SGPP::float_t>& x) {
  SGPP::float_t mean = calculateMean(x);
  SGPP::float_t var = 0.0;

  for (size_t i = 0; i < x.size(); i++) {
    var += (x[i] - mean) * (x[i] - mean);
  }

  var /= static_cast<SGPP::float_t>(x.size());
  return var;
}

void orthogonalityTest(base::DataMatrix& A) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j <= i; j++) {
      SGPP::float_t entry = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry += A(i, l) * A(j, l);
      }

      BOOST_CHECK_SMALL(entry - ((i == j) ? static_cast<SGPP::float_t>(1.0) :
                                 static_cast<SGPP::float_t>(0.0)),
                        static_cast<SGPP::float_t>(use_double_precision ? 1e-10 : 1e-6));
    }
  }
}

void symmetryTest(base::DataMatrix& A) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++) {
      BOOST_CHECK_CLOSE(A(i, j), A(j, i), 1e-10);
    }
  }
}

void similiarityTest(base::DataMatrix& A, base::DataMatrix& V,
                     base::DataMatrix& B) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());
  BOOST_CHECK_EQUAL(n, V.getNrows());
  BOOST_CHECK_EQUAL(n, V.getNcols());
  BOOST_CHECK_EQUAL(n, B.getNrows());
  BOOST_CHECK_EQUAL(n, B.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      SGPP::float_t entry1 = 0.0, entry2 = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry1 += A(i, l) * V(l, j);
        entry2 += V(i, l) * B(l, j);
      }

      BOOST_CHECK_CLOSE(entry1, entry2, (use_double_precision ? 5e-5 : 1.0));
    }
  }
}

void generateRandomMatrix(base::DataMatrix& A) {
  for (size_t i = 0; i < A.getNrows(); i++) {
    for (size_t j = 0; j < A.getNcols(); j++) {
      A(i, j) = RandomNumberGenerator::getInstance().getGaussianRN();
    }
  }
}

template <class T>
void randomMatrixEntry(T& x) {
  x = static_cast<T>(RandomNumberGenerator::getInstance().getUniformIndexRN(100));
}

template <>
void randomMatrixEntry(float& x) {
  x = static_cast<float>(RandomNumberGenerator::getInstance().getUniformRN(-100.0,
                         100.0));
}

template <>
void randomMatrixEntry(double& x) {
  x = static_cast<double>(RandomNumberGenerator::getInstance().getUniformRN(
                            -100.0, 100.0));
}

template <>
void randomMatrixEntry(std::string& x) {
  const size_t length = RandomNumberGenerator::getInstance().getUniformIndexRN(
                          100);
  x.clear();

  for (size_t i = 0; i < length; i++) {
    x += static_cast<char>(32 +
                           RandomNumberGenerator::getInstance().getUniformIndexRN(96));
  }
}

template <class T>
void testReadWriteMatrix(std::vector<T>& A1, size_t m1, size_t n1) {
  std::vector<T> A2;
  A1.clear();

  for (size_t i = 0; i < m1; i++) {
    for (size_t j = 0; j < n1; j++) {
      T entry;
      randomMatrixEntry(entry);
      A1.push_back(entry);
    }
  }

  size_t m2, n2;

  {
    const std::string fileName = "testTools_matrix.tmp";
    file_io::writeMatrix(fileName, A1, m1, n1);
    file_io::readMatrix(fileName, A2, m2, n2);
    std::remove(fileName.c_str());
  }

  BOOST_CHECK_EQUAL(m1, m2);
  BOOST_CHECK_EQUAL(n1, n2);
  BOOST_CHECK_EQUAL(A1.size(), A2.size());

  for (size_t k = 0; k < A1.size(); k++) {
    BOOST_CHECK_EQUAL(A1[k], A2[k]);
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteGrid) {
  // Test SGPP::optimization::file_io::readGrid/writeGrid.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  for (size_t d = 1; d < 6; d++) {
    std::unique_ptr<base::Grid> grid1(base::Grid::createLinearGrid(d)),
        grid2(base::Grid::createLinearGrid(d));
    std::unique_ptr<base::GridGenerator> gridGen(grid1->createGridGenerator());
    gridGen->regular(3);

    {
      const std::string fileName = "testTools_grid.tmp";
      file_io::writeGrid(fileName, *grid1->getStorage());
      file_io::readGrid(fileName, *grid2->getStorage());
      std::remove(fileName.c_str());
    }

    gridEqualityTest(*grid1, *grid2);

    base::DataVector functionValues1(grid1->getStorage()->size());
    base::DataVector functionValues2(0);

    for (size_t k = 0; k < functionValues1.getSize(); k++) {
      functionValues1[k] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    {
      const std::string fileName = "testTools_grid.tmp";
      file_io::writeGrid(fileName, *grid1->getStorage(), functionValues1);
      std::unique_ptr<base::Grid> grid2(base::Grid::createLinearGrid(d));
      file_io::readGrid(fileName, *grid2->getStorage(), functionValues2);
      std::remove(fileName.c_str());
    }

    gridEqualityTest(*grid1, *grid2);
    BOOST_CHECK_EQUAL(functionValues1.getSize(), functionValues2.getSize());

    for (size_t k = 0; k < functionValues1.getSize(); k++) {
      BOOST_CHECK_EQUAL(functionValues1[k], functionValues2[k]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteMatrix) {
  // Test SGPP::optimization::file_io::readMatrix/writeMatrix.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t m1 = 100;
  const size_t n1 = 200;

  // test read/write with std::vector<T>
  {
    std::vector<SGPP::float_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<float> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<double> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint8_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint16_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint32_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint64_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<std::string> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  // test read/write with DataMatrix
  {
    base::DataMatrix A1(m1, n1), A2(0, 0);

    for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
        A1(i, j) = RandomNumberGenerator::getInstance().getUniformRN();
      }
    }

    {
      const std::string fileName = "testTools_matrix.tmp";
      file_io::writeMatrix(fileName, A1);
      file_io::readMatrix(fileName, A2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(A1.getNrows(), A2.getNrows());
    BOOST_CHECK_EQUAL(A1.getNcols(), A2.getNcols());

    for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
        BOOST_CHECK_EQUAL(A1(i, j), A2(i, j));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteVector) {
  // Test SGPP::optimization::file_io::readVector/writeVector.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 100;

  // test read/write with std::vector<SGPP::float_t>
  {
    std::vector<SGPP::float_t> v1, v2;

    for (size_t i = 0; i < n; i++) {
      v1.push_back(RandomNumberGenerator::getInstance().getUniformRN());
    }

    {
      const std::string fileName = "testTools_matrix.tmp";
      file_io::writeVector(fileName, v1);
      file_io::readVector(fileName, v2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(v1.size(), v2.size());

    for (size_t k = 0; k < v1.size(); k++) {
      BOOST_CHECK_EQUAL(v1[k], v2[k]);
    }
  }

  // test read/write with DataVector
  {
    base::DataVector v1(n), v2(0);

    for (size_t i = 0; i < n; i++) {
      v1[i] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    {
      const std::string fileName = "testTools_vector.tmp";
      file_io::writeVector(fileName, v1);
      file_io::readVector(fileName, v2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(v1.getSize(), v2.getSize());

    for (size_t i = 0; i < v1.getSize(); i++) {
      BOOST_CHECK_EQUAL(v1[i], v2[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestRandomNumberGenerator) {
  // Test SGPP::optimization::RandomNumberGenerator.
  const size_t seed = 42;
  const size_t N = 20000;
  std::vector<SGPP::float_t> numbers(N);

  // set and test seed getting/setting
  RandomNumberGenerator::getInstance().setSeed();
  RandomNumberGenerator::getInstance().setSeed(seed);
  BOOST_CHECK_EQUAL(RandomNumberGenerator::getInstance().getSeed(), seed);

  // test continuous uniform random numbers
  {
    for (size_t i = 0; i < N; i++) {
      numbers[i] = RandomNumberGenerator::getInstance().getUniformRN();
      BOOST_CHECK_GE(numbers[i], 0.0);
      BOOST_CHECK_LE(numbers[i], 1.0);
    }

    BOOST_CHECK_SMALL(calculateMean(numbers) - static_cast<SGPP::float_t>(0.5),
                      static_cast<SGPP::float_t>(
                        use_double_precision ? 1e-3 : 1e-2));
    BOOST_CHECK_SMALL(calculateVariance(numbers) -
                      static_cast<SGPP::float_t>(1.0) /
                      static_cast<SGPP::float_t>(12.0),
                      static_cast<SGPP::float_t>(1e-3));
  }

  // test Gaussian random numbers
  {
    std::vector<SGPP::float_t> mus = {0.0, 12.3, -42.0, 13.37};
    std::vector<SGPP::float_t> sigmas = {1.0, 2.6, 8.1, 0.3};

    for (size_t k = 0; k < mus.size(); k++) {
      for (size_t i = 0; i < N; i++) {
        numbers[i] = RandomNumberGenerator::getInstance().getGaussianRN(sigmas[k],
                     mus[k]);
      }

      BOOST_CHECK_SMALL(calculateMean(numbers) - mus[k],
                        static_cast<SGPP::float_t>(0.1) * sigmas[k]);
      BOOST_CHECK_SMALL(calculateVariance(numbers) - sigmas[k] * sigmas[k],
                        static_cast<SGPP::float_t>(0.1) *
                        sigmas[k] * sigmas[k]);
    }
  }

  // test discrete uniform random numbers
  for (size_t k = 1; k < 11; k++) {
    for (size_t i = 0; i < N; i++) {
      numbers[i] = static_cast<SGPP::float_t>(
                     RandomNumberGenerator::getInstance().getUniformIndexRN(k));
      BOOST_CHECK_EQUAL(numbers[i], static_cast<int>(numbers[i]));
      BOOST_CHECK_GE(numbers[i], 0);
      BOOST_CHECK_LE(numbers[i], k - 1);
    }

    SGPP::float_t kDbl = static_cast<SGPP::float_t>(k);
    BOOST_CHECK_SMALL(calculateMean(numbers) -
                      (kDbl - static_cast<SGPP::float_t>(1.0)) /
                      static_cast<SGPP::float_t>(2.0),
                      static_cast<SGPP::float_t>(0.01) * kDbl);
    BOOST_CHECK_SMALL(calculateVariance(numbers) -
                      (kDbl * kDbl - static_cast<SGPP::float_t>(1.0)) /
                      static_cast<SGPP::float_t>(12.0),
                      static_cast<SGPP::float_t>(0.01) * kDbl * kDbl);
  }
}

BOOST_AUTO_TEST_CASE(TestHouseholderTransformation) {
  // Test SGPP::optimization::math::householderTransformation.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  base::DataMatrix A(n, n);
  generateRandomMatrix(A);

  const std::vector<size_t> ps = {0, 2, 5};
  const std::vector<size_t> qs = {0, 1, 9};

  for (size_t k = 0; k < ps.size(); k++) {
    const size_t p = ps[k];
    const size_t q = qs[k];
    const size_t m = n - p;
    base::DataMatrix Q(m, m);
    math::householderTransformation(A, p, q, Q);
    symmetryTest(Q);
    orthogonalityTest(Q);

    for (size_t i = p + 1; i < n; i++) {
      SGPP::float_t entry = 0.0;

      for (size_t l = p; l < n; l++) {
        entry += Q(i - p, l - p) * A(l, q);
      }

      BOOST_CHECK_SMALL(entry, static_cast<SGPP::float_t>(
                          use_double_precision ? 1e-10 : 1e-6));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestHessenbergForm) {
  // Test SGPP::optimization::math::hessenbergForm.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  base::DataMatrix A(n, n);
  generateRandomMatrix(A);
  base::DataMatrix H(A);
  base::DataMatrix V(n, n);
  math::hessenbergForm(H, V);

  orthogonalityTest(V);
  similiarityTest(A, V, H);

  for (size_t i = 2; i < n; i++) {
    for (size_t j = 0; j < i - 1; j++) {
      BOOST_CHECK_SMALL(H(i, j), static_cast<SGPP::float_t>(1e-10));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestQRDecomposition) {
  // Test SGPP::optimization::math::QRDecomposition.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  base::DataMatrix A(n, n);
  generateRandomMatrix(A);
  base::DataMatrix R(A);
  base::DataMatrix Q(n, n);
  math::QRDecomposition(R, Q);

  orthogonalityTest(Q);

  for (size_t i = 1; i < n; i++) {
    for (size_t j = 0; j < i - 1; j++) {
      BOOST_CHECK_SMALL(R(i, j), static_cast<SGPP::float_t>(1e-10));
    }
  }

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      SGPP::float_t entry = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry += Q(i, l) * R(l, j);
      }

      BOOST_CHECK_CLOSE(A(i, j), entry, static_cast<SGPP::float_t>(
                          use_double_precision ? 1e-10 : 1e-2));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestSchurDecomposition) {
  // Test SGPP::optimization::math::schurDecomposition.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  base::DataMatrix A(n, n);
  generateRandomMatrix(A);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++) {
      A(i, j) = A(j, i);
    }
  }

  base::DataMatrix S(A);
  base::DataMatrix V(n, n);
  math::schurDecomposition(S, V);

  orthogonalityTest(V);
  similiarityTest(A, V, S);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
        BOOST_CHECK_SMALL(S(i, j), static_cast<SGPP::float_t>(
                            use_double_precision ? 1e-8 : 1e-5));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestPrinter) {
  // Test SGPP::optimization::Printer.

  // redirect std::cout to not confuse the user
  const std::string fileName = "testTools_Printer::getInstance().tmp";
  std::ofstream outStream(fileName);
  Printer::getInstance().setStream(&outStream);

  Printer::getInstance().setVerbosity(2);
  Printer::getInstance().printStatusBegin("Testing status printing...");
  Printer::getInstance().printStatusUpdate("Test status update 1");
  Printer::getInstance().printStatusUpdate("Test status update 2");
  Printer::getInstance().printStatusEnd("Testing status printing ended.");

  Printer::getInstance().getMutex();

  const SGPP::float_t duration = Printer::getInstance().getLastDurationSecs();
  BOOST_CHECK_GE(duration, 0.0);
  BOOST_CHECK_LE(duration, 0.01);

  base::DataMatrix A(3, 3, 0.0);
  A(0, 1) = 12.3;
  A(1, 2) = 42.1337;
  FullSLE sle(A);
  Printer::getInstance().printSLE(sle);

  const size_t d = 1;
  const size_t p = 5;
  const size_t N = 10;
  test_problems::SphereObjective f(d);
  std::unique_ptr<base::Grid> grid(base::Grid::createModBsplineGrid(d, p));
  IterativeGridGeneratorRitterNovak gridGen(f, *grid, N, 0.85);
  BOOST_CHECK(gridGen.generate());
  Printer::getInstance().printIterativeGridGenerator(gridGen);

  // undo redirection
  outStream.close();
  std::remove(fileName.c_str());
  Printer::getInstance().setStream(&std::cout);
}
