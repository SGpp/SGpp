#define BOOST_TEST_DYN_LINK
#include <vector>
#include <boost/test/unit_test.hpp>

#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

using namespace SGPP::base;

BOOST_AUTO_TEST_SUITE(TestGridFactory)

BOOST_AUTO_TEST_CASE(testCreation) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridStorage* storage = factory->getStorage();
  BOOST_CHECK( storage != NULL );


  delete( factory );
}


BOOST_AUTO_TEST_CASE(testSerializationLinear) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationModLinear) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createModLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundary) {
  // Uses Linear grid for tests
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  delete( factory );
  delete( gen );
  delete( newfac );
}



BOOST_AUTO_TEST_CASE(testSerializationLinearBoundary) {
  // Uses Linear grid for tests
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationPrewavelet) {
  // Uses Linear grid for tests
  Grid* factory = Grid::createPrewaveletGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearBoundingBox) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  {
    BoundingBox* boundingBox = factory->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox->setBoundary( 0, tempBound );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  {
    BoundingBox* boundingBox = newfac->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    BOOST_CHECK( 0.0 == tempBound.leftBoundary );
    BOOST_CHECK( 100.0 == tempBound.rightBoundary );
    BOOST_CHECK( false == tempBound.bDirichletLeft );
    BOOST_CHECK( false == tempBound.bDirichletRight );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationModLinearBoundingBox) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createModLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  {
    BoundingBox* boundingBox = factory->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox->setBoundary( 0, tempBound );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  {
    BoundingBox* boundingBox = newfac->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    BOOST_CHECK( 0.0 == tempBound.leftBoundary );
    BOOST_CHECK( 100.0 == tempBound.rightBoundary );
    BOOST_CHECK( false == tempBound.bDirichletLeft );
    BOOST_CHECK( false == tempBound.bDirichletRight );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundaryBoundingBox) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  {
    BoundingBox* boundingBox = factory->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox->setBoundary( 0, tempBound );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  {
    BoundingBox* boundingBox = newfac->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    BOOST_CHECK( 0.0 == tempBound.leftBoundary );
    BOOST_CHECK( 100.0 == tempBound.rightBoundary );
    BOOST_CHECK( false == tempBound.bDirichletLeft );
    BOOST_CHECK( false == tempBound.bDirichletRight );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearBoundaryBoundingBox) {
  // Uses Linear grid for tests

  Grid* factory = Grid::createLinearBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  {
    BoundingBox* boundingBox = factory->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    tempBound.leftBoundary = 0.0;
    tempBound.rightBoundary = 100.0;
    tempBound.bDirichletLeft = false;
    tempBound.bDirichletRight = false;
    boundingBox->setBoundary( 0, tempBound );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  {
    BoundingBox* boundingBox = newfac->getBoundingBox();
    DimensionBoundary tempBound = boundingBox->getBoundary(0);
    BOOST_CHECK( 0.0 == tempBound.leftBoundary );
    BOOST_CHECK( 100.0 == tempBound.rightBoundary );
    BOOST_CHECK( false == tempBound.bDirichletLeft );
    BOOST_CHECK( false == tempBound.bDirichletRight );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  Grid* factory = Grid::createLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);


  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    srcLeaf.push_back( factory->getStorage()->get(i)->isLeaf() );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    BOOST_CHECK( newfac->getStorage()->get(i)->isLeaf() == srcLeaf[i] );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationModLinearWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  Grid* factory = Grid::createModLinearGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    srcLeaf.push_back( factory->getStorage()->get(i)->isLeaf() );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    BOOST_CHECK( newfac->getStorage()->get(i)->isLeaf() == srcLeaf[i] );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearTruncatedBoundaryWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    srcLeaf.push_back( factory->getStorage()->get(i)->isLeaf() );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    BOOST_CHECK( newfac->getStorage()->get(i)->isLeaf() == srcLeaf[i] );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}


BOOST_AUTO_TEST_CASE(testSerializationLinearBoundaryWithLeaf) {
  // Uses Linear grid for tests

  std::vector<bool> srcLeaf;

  Grid* factory = Grid::createLinearBoundaryGrid(2);
  BOOST_CHECK( factory != NULL );

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(3);

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    srcLeaf.push_back( factory->getStorage()->get(i)->isLeaf() );
  }

  std::string str = factory->serialize();
  BOOST_CHECK( str.size() > 0 );

  Grid* newfac = Grid::unserialize( str );
  BOOST_CHECK( newfac != NULL );

  BOOST_CHECK( factory->getStorage()->size() == newfac->getStorage()->size() );

  for ( size_t i = 0; i < factory->getStorage()->size(); ++i ) {
    BOOST_CHECK( newfac->getStorage()->get(i)->isLeaf() == srcLeaf[i] );
  }

  delete( factory );
  delete( gen );
  delete( newfac );
}

// end test suite TestGridFactory
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(TestLinearGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  Grid* factory = Grid::createLinearGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  BOOST_CHECK( gen != NULL );

  BOOST_CHECK( storage->size() == 0 );
  gen->regular(3);
  BOOST_CHECK( storage->size() == 17 );

  // This should fail
  BOOST_CHECK_THROW( gen->regular(3), generation_exception );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement) {
  Grid* factory = Grid::createLinearGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  BOOST_CHECK( storage->size() == 1 );
  DataVector alpha(1);
  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 5 );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  Grid* factory = Grid::createLinearGrid(1);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(2);

  DataVector alpha( factory->getStorage()->size() );
  DataMatrix p(1, 1);
  DataVector beta(1);

  alpha.setAll(0.0);
  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  OperationMultipleEval* opb = SGPP::op_factory::createOperationMultipleEval( *factory, p );
  opb->multTranspose(beta, alpha);

  BOOST_CHECK_CLOSE( alpha[0], 0.5, 0.0 );
  BOOST_CHECK_CLOSE( alpha[1], 1.0, 0.0 );
  BOOST_CHECK_CLOSE( alpha[2], 0.0, 0.0 );

  alpha.setAll(0.0);
  alpha[0] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult(alpha, beta);
  BOOST_CHECK_CLOSE( beta[0], 0.5, 0.0 );

  delete( gen );
  delete( opb );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  Grid* factory = Grid::createLinearGrid(1);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  DataVector alpha( factory->getStorage()->size() );
  alpha.setAll(1.0);

  DataVector p(1);
  p.setAll(0.25);

  OperationEval* eval = SGPP::op_factory::createOperationEval( *factory );
  BOOST_CHECK_CLOSE( eval->eval(alpha, p), 0.5, 0.0 );

  delete( gen );
  delete( eval );
  delete( factory );
}

// end test suite TestLinearGrid
BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(TestLinearTruncatedBoundaryGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  BOOST_CHECK( gen != NULL );

  BOOST_CHECK( storage->size() == 0 );
  gen->regular(3);
  BOOST_CHECK( storage->size() == 49 );

  // This should fail
  BOOST_CHECK_THROW( gen->regular(3), generation_exception );

  delete factory;
  delete gen;
}

BOOST_AUTO_TEST_CASE(testRefinement2d) {
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  BOOST_CHECK( storage->size() == 9 );
  DataVector alpha( 9 );
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  alpha[3] = 0.0;
  alpha[4] = 0.0;
  alpha[5] = 0.0;
  alpha[6] = 0.0;
  alpha[7] = 0.0;
  alpha[8] = 1.0;

  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 21 );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(3);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  BOOST_CHECK( storage->size() == 27 );
  DataVector alpha( 27 );
  alpha.setAll(0.0);

  alpha[26] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 81 );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement3d) {
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(1);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(2);

  DataVector alpha( factory->getStorage()->size() );
  DataMatrix p(1, 1);
  DataVector beta(1);

  alpha.setAll(0.0);
  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  OperationMultipleEval* opb = SGPP::op_factory::createOperationMultipleEval( *factory, p);
  opb->multTranspose(beta, alpha);

  BOOST_CHECK_CLOSE( alpha[0], 0.75, 0.0 );
  BOOST_CHECK_CLOSE(alpha[1], 0.25, 0.0 );
  BOOST_CHECK_CLOSE(alpha[2], 0.5, 0.0 );
  BOOST_CHECK_CLOSE(alpha[3], 1.0, 0.0 );
  BOOST_CHECK_CLOSE(alpha[4], 0.0, 0.0 );

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult(alpha, beta);
  BOOST_CHECK_CLOSE(beta[0], 0.5, 0.0 );

  delete( gen );
  delete( opb );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  Grid* factory = Grid::createLinearTruncatedBoundaryGrid(1);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  DataVector alpha( factory->getStorage()->size() );
  alpha.setAll(1.0);

  DataVector p(1);
  p.setAll(0.25);

  OperationEval* eval = SGPP::op_factory::createOperationEval( *factory );

  BOOST_CHECK_CLOSE( eval->eval(alpha, p), 1.5, 0.0 );

  delete( gen );
  delete( eval );
  delete( factory );
}

// end test suite TestLinearTruncatedBoundaryGrid
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(TestLinearBoundaryGrid)

BOOST_AUTO_TEST_CASE(testGeneration) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  BOOST_CHECK( gen != NULL );

  BOOST_CHECK( storage->size() == 0 );
  gen->regular(3);
  BOOST_CHECK( storage->size() == 37 );

  // This should fail
  BOOST_CHECK_THROW( gen->regular(3), generation_exception );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement2d_one) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(0);

  BOOST_CHECK( storage->size() == 4);

  DataVector alpha(4);
  alpha.setAll(0.0);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 8);


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement2d_two) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(0);

  DataVector alpha(4);
  alpha.setAll(0.0);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );

  DataVector alpha2( 8 );
  alpha2.setAll( 0.0 );

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2( &alpha2 );

  gen->refine( &func2 );
  BOOST_CHECK( storage->size() == 13 );


  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement2d_three) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(0);

  DataVector alpha(4);
  alpha.setAll(0.0);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );

  DataVector alpha2( 8 );
  alpha2.setAll( 0.0 );

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2( &alpha2 );

  gen->refine( &func2 );

  DataVector alpha3(13);
  alpha3.setAll(0.0);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3( &alpha3 );

  gen->refine( &func3 );
  BOOST_CHECK( storage->size() == 18 );

  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement2d_four) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(0);

  DataVector alpha(4);
  alpha.setAll(0.0);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );

  DataVector alpha2( 8 );
  alpha2.setAll( 0.0 );

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2( &alpha2 );

  gen->refine( &func2 );

  DataVector alpha3(13);
  alpha3.setAll(0.0);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3( &alpha3 );

  gen->refine( &func3 );

  DataVector alpha4( 18 );
  alpha4.setAll(0.0);

  alpha4[12] = 1.0;
  SurplusRefinementFunctor func4( &alpha4 );

  gen->refine( &func4 );
  BOOST_CHECK( storage->size() == 25 );

  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testRefinement2d_five) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(0);

  DataVector alpha(4);
  alpha.setAll(0.0);

  alpha[0] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );

  DataVector alpha2( 8 );
  alpha2.setAll( 0.0 );

  alpha2[4] = 1.0;
  SurplusRefinementFunctor func2( &alpha2 );

  gen->refine( &func2 );

  DataVector alpha3(13);
  alpha3.setAll(0.0);

  alpha3[11] = 1.0;
  SurplusRefinementFunctor func3( &alpha3 );

  gen->refine( &func3 );

  DataVector alpha4( 18 );
  alpha4.setAll(0.0);

  alpha4[12] = 1.0;
  SurplusRefinementFunctor func4( &alpha4 );

  gen->refine( &func4 );

  DataVector alpha5( 25 );
  alpha5.setAll(0.0);

  alpha5[23] = 1.0;
  SurplusRefinementFunctor func5( &alpha5 );

  gen->refine( &func5 );

  BOOST_CHECK( storage->size() == 29 );

  delete( gen );
  delete( factory );
}


BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(2);

  DataVector alpha( factory->getStorage()->size() );
  DataMatrix p(1, 1);
  DataVector beta(1);

  alpha.setAll(0.0);
  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  OperationMultipleEval* opb = SGPP::op_factory::createOperationMultipleEval( *factory, p );
  opb->multTranspose( beta, alpha );

  BOOST_CHECK_CLOSE( alpha[0], 0.75, 0.0 );
  BOOST_CHECK_CLOSE( alpha[1], 0.25, 0.0 );
  BOOST_CHECK_CLOSE( alpha[2], 0.5, 0.0 );
  BOOST_CHECK_CLOSE( alpha[3], 1.0, 0.0 );
  BOOST_CHECK_CLOSE( alpha[4], 0.0, 0.0 );

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  p.set(0, 0, 0.25);

  beta[0] = 0.0;

  opb->mult( alpha, beta);
  BOOST_CHECK_CLOSE( beta[0], 0.5, 0.0 );

  delete gen;
  delete opb;
  delete factory;
}


BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  Grid* factory = Grid::createLinearBoundaryGrid(2);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  DataVector alpha( factory->getStorage()->size() );
  alpha.setAll( 1.0 );

  DataVector p( 1 );
  p.setAll(0.25);

  OperationEval* eval = SGPP::op_factory::createOperationEval( *factory );

  BOOST_CHECK_CLOSE( eval->eval( alpha, p ), 1.5, 0.0 );

  delete( gen );
  delete( eval );
  delete( factory );
}

// end test suite TestLinearBoundaryGrid
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( TestLinearStretchedTruncatedBoundaryGrid )

BOOST_AUTO_TEST_CASE(testGeneration) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;

  DimensionBoundary dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;

  std::vector<DimensionBoundary> dimBoundaryVector(2);
  dimBoundaryVector[0] = dimBound;
  dimBoundaryVector[1] = dimBound;

  std::vector<Stretching1D> str1dvector(2);
  str1dvector[0] = str1d;
  str1dvector[1] = str1d;
  Stretching stretch(2, dimBoundaryVector, str1dvector);

  Grid* factory = Grid::createLinearStretchedTruncatedBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  BOOST_CHECK( gen != NULL );

  BOOST_CHECK( storage->size() == 0 );
  gen->regular(3);
  BOOST_CHECK( storage->size() == 49 );

  // This should fail
  BOOST_CHECK_THROW( gen->regular(3), generation_exception );

  delete gen;
  delete factory;
}


BOOST_AUTO_TEST_CASE(testRefinement2d) {
  Grid* factory = Grid::createLinearStretchedTruncatedBoundaryGrid(2);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  BOOST_CHECK( storage->size() == 9 );
  DataVector alpha(9);
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  alpha[3] = 0.0;
  alpha[4] = 0.0;
  alpha[5] = 0.0;
  alpha[6] = 0.0;
  alpha[7] = 0.0;
  alpha[8] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 21 );

  delete gen;
  delete factory;
}


BOOST_AUTO_TEST_CASE(testRefinement3d) {
  Grid* factory = Grid::createLinearStretchedTruncatedBoundaryGrid(3);
  GridStorage* storage = factory->getStorage();

  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  BOOST_CHECK( storage->size() == 27 );
  DataVector alpha(27);
  alpha.setAll(0.0);

  alpha[26] = 1.0;
  SurplusRefinementFunctor func( &alpha );

  gen->refine( &func );
  BOOST_CHECK( storage->size() == 81 );

  delete gen;
  delete factory;
}


BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;
  DimensionBoundary dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch( 1, &dimBound, &str1d );

  Grid* factory = Grid::createLinearStretchedTruncatedBoundaryGrid(1);
  factory->getStorage()->setStretching(stretch);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(2);

  DataVector alpha( factory->getStorage()->size() );

  DataMatrix p(1, 1);
  DataVector beta(1);

  alpha.setAll(0.0);
  p.set(0, 0, 0.25);
  beta[0] = 1.0;

  OperationMultipleEval* opb = SGPP::op_factory::createOperationMultipleEval( *factory, p );
  opb->multTranspose( beta, alpha );

#if USE_DOUBLE_PRECISION == 1
  BOOST_CHECK_CLOSE(alpha[0], SGPP::float_t(1.038461538), SGPP::float_t(1e-6) );
  BOOST_CHECK_CLOSE(alpha[1], SGPP::float_t(-0.038461538461538464), SGPP::float_t(1e-6) );
  BOOST_CHECK_CLOSE(alpha[2], SGPP::float_t(-0.18237143795284394), SGPP::float_t(1e-6) );
  BOOST_CHECK_CLOSE(alpha[3], SGPP::float_t(-0.53513915), SGPP::float_t(1e-6) );
  BOOST_CHECK_CLOSE(alpha[4], SGPP::float_t(0.0), SGPP::float_t(1e-6) );
#else
  BOOST_CHECK_CLOSE(alpha[0], SGPP::float_t(1.038461538), SGPP::float_t(1e-5) );
  BOOST_CHECK_CLOSE(alpha[1], SGPP::float_t(-0.038461538461538464), SGPP::float_t(1e-5) );
  BOOST_CHECK_CLOSE(alpha[2], SGPP::float_t(-0.18237143795284394), SGPP::float_t(1e-4) );
  BOOST_CHECK_CLOSE(alpha[3], SGPP::float_t(-0.53513915), SGPP::float_t(1e-4) );
  BOOST_CHECK_CLOSE(alpha[4], SGPP::float_t(0.0), SGPP::float_t(1e-6) );
#endif

  alpha.setAll(0.0);
  alpha[2] = 1.0;

  p.set( 0, 0, 0.25 );

  beta[0] = 0.0;

  opb->mult( alpha, beta );
#if USE_DOUBLE_PRECISION == 1
  BOOST_CHECK_CLOSE( beta[0], SGPP::float_t(-0.182371437), SGPP::float_t(1e-6) );
#else
  BOOST_CHECK_CLOSE( beta[0], SGPP::float_t(-0.182371437), SGPP::float_t(1e-4) );
#endif

  delete gen;
  delete factory;
}


BOOST_AUTO_TEST_CASE(testOperationEval_eval) {
  /*
  from pysgpp import Grid, DataVector, Stretching, Stretching1D, DimensionBoundary

  str1d = Stretching1D()
  str1d.type='log'
  str1d.x_0=1
  str1d.xsi=10
  dimBound = DimensionBoundary()
  dimBound.leftBoundary=0.5
  dimBound.rightBoundary=7
  stretch=Stretching(1,dimBound,str1d)

  factory = Grid.createLinearStretchedTruncatedBoundaryGrid(1)
  factory.getStorage().setStretching(stretch)
  gen = factory.createGridGenerator()
  gen.regular(1)

  alpha = DataVector(factory.getStorage().size())
  alpha.setAll(1.0)

  p = DataVector(1)
  p.setAll(0.25)

  eval = createOperationEval(factory)

  self.failUnlessAlmostEqual(eval.eval(alpha, p), 0.8176285620)
  */
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;
  DimensionBoundary dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch( 1, &dimBound, &str1d );

  Grid* factory = Grid::createLinearStretchedTruncatedBoundaryGrid(1);
  factory->getStorage()->setStretching(stretch);
  GridGenerator* gen = factory->createGridGenerator();
  gen->regular(1);

  DataVector alpha( factory->getStorage()->size() );
  alpha.setAll(1.0);

  DataVector p(1);
  p.setAll(0.25);

  OperationEval* eval = SGPP::op_factory::createOperationEval( *factory );

#if USE_DOUBLE_PRECISION == 1
  BOOST_CHECK_CLOSE( eval->eval( alpha, p ), SGPP::float_t(0.8176285620), SGPP::float_t(1e-8) );
#else
  BOOST_CHECK_CLOSE( eval->eval( alpha, p ), SGPP::float_t(0.8176285620), SGPP::float_t(1e-5) );
#endif

  delete gen;
  delete factory;
}
// end test suite TestLinearStretchedTruncatedBoundaryGrid
BOOST_AUTO_TEST_SUITE_END()



