// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <stdlib.h>

#include <sgpp/combigrid/multigridFG/utils/RunTikhonov.hpp>
#include <sgpp/combigrid/multigridFG/operators/TikhonovOperator.hpp>
#include <sgpp/combigrid/multigridFG/multigrid/Multigrid.hpp>

using namespace combigrid;

void combigrid::RunTikhonov::readInInput(
  const std::string& XfileS ,
  const std::string& YFileS ,
  int& dimensions ,  int& nrPoints ,
  std::vector<double>& XCoords,
  std::vector<double>& YPoint) {

  int verb = 3;
  dimensions = 0;
  nrPoints = 0;

  std::ifstream Xfile;
  std::ifstream Yfile;

  Xfile.open( XfileS.c_str() , std::ios::in );
  Yfile.open( YFileS.c_str() , std::ios::in );

  // the string which we read in one line
  std::string str_tmp;
  double d1;
  int i;
  char* pEnd, *pStart;

  // read in one line to detect the dimension
  getline( Xfile , str_tmp );
  dimensions = 0;
  pStart = &(str_tmp[0]); //str_tmp.c_str();
  std::cout << " first line: " << str_tmp << std::endl;

  while ( pStart[0] != 0 ) {
    d1 = strtod (pStart , &pEnd);
    std::cout << " number: " << d1 << " pEnd:" << (int)(pEnd[0]) << std::endl;
    pStart = pEnd;

    if (pStart != NULL) {
      dimensions++;
      XCoords.push_back(d1);
    }
  }

  std::cout << " dim: " << dimensions << std::endl;

  // now we have the dimension, read in the X
  getline( Xfile , str_tmp );

  while (!Xfile.eof()) {
    // for each line get the
    pStart = &(str_tmp[0]); //str_tmp.c_str();

    for (i = 0 ; i < dimensions ; i++) {
      d1 = strtod (pStart , &pEnd);
      pStart = pEnd;
      XCoords.push_back(d1);
    }

    getline( Xfile , str_tmp );
  }

  // read in Y
  getline( Yfile , str_tmp );

  while (!Yfile.eof()) {
    // there is only one number per line
    pStart = &(str_tmp[0]); //str_tmp.c_str();
    d1 = strtod (pStart , &pEnd);
    YPoint.push_back(d1);
    getline( Yfile , str_tmp );
    nrPoints++;
  }

  // print the results
  int j = 0;

  if (verb > 3) {
    for ( i = 0 ; i < nrPoints ; i++) {
      std::cout << " X[" << i << "] ( " << XCoords[i * dimensions];

      for ( j = 1 ; j < dimensions ; j++) {
        std::cout << "," << XCoords[i * dimensions + j];
      }

      std::cout << " ) Y = " << YPoint[i] << std::endl;
    }
  }
}

FullGridD* combigrid::RunTikhonov::computeFGTikhonov(
  combigrid::GridDomain& domain ,
  const std::vector<int>& levels,
  double lambda,
  const std::string& Xfile ,
  const std::string& YFile ) {

  int dimensions = 0;
  int nrPoints = 0;
  int verb = 6;
  std::vector<double> XCoords;
  std::vector<double> YPoint;
  std::vector<double> unknowns;
  FullGridD* fg;
  // first read in the results
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... read INPUT");
  combigrid::RunTikhonov::readInInput(Xfile , YFile , dimensions , nrPoints , XCoords, YPoint);

  // create the full grid
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create Full Grid");
  fg = new FullGridD( dimensions , levels );
  fg->createFullGrid();
  fg->setDomain( &domain );

  // create the Tikhonov operator
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::TikhonovOperator");
  combigrid::TikhonovOperator op( fg , nrPoints, lambda , &XCoords , &YPoint );

  // create the multigrid method
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::Multigrid");
  combigrid::Multigrid multigrid( &(op) , fg );

  // solve using only the smoother
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... solve multigird");
  unknowns.resize(fg->getNrElements(), 0.0);
  multigrid.solveCS( unknowns , 1e-8 , false );
  //multigrid.solveSmoothing( unknowns , 1e-8);
  //multigrid.solveCG( unknowns , 1e-6 );
  //multigridFAS.solveFAS( unknowns , 1e-8 );

  // copy the solution back
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... write solution back and return");

  for (int i = 0 ; i < fg->getNrElements() ; i++) {
    fg->getElementVector()[i] = unknowns[i];
    //COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
  }

  // return full grid
  return fg;
}

void combigrid::RunTikhonov::computeFGTikhonov_FG(
  FullGridD* fg,
  std::vector<double>& unknowns ,
  double lambda,
  std::vector<double>& XCoords,
  std::vector<double>& YPoint ) {

  //int dimensions = (int) ::round( ((double)XCoords.size()) / ((double)YPoint.size()) );
  int nrPoints = static_cast<int>( YPoint.size() );
  int verb = 6;

  // create the Tikhonov operator
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::TikhonovOperator");
  combigrid::TikhonovOperator op( fg , nrPoints, lambda , &XCoords , &YPoint );

  // create the multigrid method
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... create combigrid::Multigrid");
  combigrid::Multigrid multigrid( &(op) , fg );

  // solve using only the smoother
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... solve multigird");
  multigrid.solveCS( unknowns , 1e-8 , false );
  //multigrid.solveSmoothing( unknowns , 1e-8);
  //multigrid.solveCG( unknowns , 1e-6 );
  //multigridFAS.solveFAS( unknowns , 1e-8 );

  // copy the solution back
  COMBIGRID_OUT_LEVEL2( verb , "RunTikhonov::computeFGTikhonov ... write solution back and return");

  for (int i = 0 ; i < fg->getNrElements() ; i++) {
    fg->getElementVector()[i] = unknowns[i];
    //COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
  }
}


void combigrid::RunTikhonov::computeTikhonov_CG(
  AbstractCombiGrid* combiG,
  double lambda,
  std::vector<double>& XCoords,
  std::vector<double>& YPoint) {

  // for each full grid call the
  for (int gr = 0 ; gr < combiG->getNrFullGrid() ; gr++) {
    computeFGTikhonov_FG( combiG->getFullGrid(gr) ,
                          combiG->getFullGrid(gr)->getElementVector() , lambda, XCoords, YPoint );
  }
}


void combigrid::RunTikhonov::computeTikhonov_FG_crossvalidation(
  FullGridD* fg,
  double& lambda,
  std::vector<double>& XCoords,
  std::vector<double>& YPoint) {

  int verb = 6;
  int trainPointNr = 3 * static_cast<int>( YPoint.size() ) / 4;
  int testPointNr = static_cast<int>( YPoint.size() ) - trainPointNr;
  int dim = (int)::round((double)XCoords.size() / (double)YPoint.size());

  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_FG_crossvalidation    , trainPointNr="
                        << trainPointNr << " , testPointNr=" << testPointNr);

  std::vector<double> XCoordsTrain(trainPointNr * dim);
  std::vector<double> YPointTrain(trainPointNr);
  std::vector<double> XCoordsTest(testPointNr * dim);
  std::vector<double> YPointTest(testPointNr);

  // copy the test and train values in special arrays
  for (int mp = 0 ; mp < trainPointNr * dim ; mp ++) {
    XCoordsTrain[mp] = XCoords[mp];
  }

  for (int mp = trainPointNr * dim ; mp < (testPointNr + trainPointNr)*dim ; mp ++) {
    XCoordsTest[mp - trainPointNr * dim] = XCoords[mp];
  }

  for (int mp = 0 ; mp < trainPointNr ; mp ++) {
    YPointTrain[mp] = YPoint[mp];
  }

  for (int mp = trainPointNr ; mp < (testPointNr + trainPointNr) ; mp ++) {
    YPointTest[mp - trainPointNr] = YPoint[mp];
  }

  double lamdba1 , lamdba2 = 2e-3 , err1 , err2 , step = (0.9) * lamdba2;
  int direction = 1 , iterNr = 0;
  lamdba1 = lamdba2;
  // we compute the initial error
  computeFGTikhonov_FG( fg, fg->getElementVector() , lamdba1 , XCoordsTrain, YPointTrain);
  err1 = measureErrorFG( fg, dim  , XCoordsTest, YPointTest );
  err2 = err1;
  lamdba2 = lamdba1;
  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_FG_crossvalidation    , initial error="
                        << err1 << " , lamdba=" << lamdba1 << " , step = " << step);

  // iterate till we find the solution , this is a simple minimum search algorithm
  // we either look in one direction or into the other, and we change the "step" size so we should hit the minimum
  while ( (iterNr < 15) && (step > 1e-8) ) {
    //
    lamdba1 = lamdba1 - ((double)direction) * step;
    computeFGTikhonov_FG( fg, fg->getElementVector() , lamdba1 , XCoordsTrain, YPointTrain);
    err1 = measureErrorFG( fg, dim  , XCoordsTest, YPointTest );
    COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_FG_crossvalidation   try step lambda=" << lamdba1 <<
                          "  error=" << err1 << " , step=" << step);

    // the step was successful
    if (err1 < err2) {
      err2 = err1;
      lamdba2 = lamdba1;

      if ( direction == 1) {
        step = step / 10.0;
      }
    } else {
      // the step leads to a higher value, change direction
      err2 = err1;
      lamdba2 = lamdba1;

      if (direction == 1) {
        step = 3.0 * step;
      }

      if (direction < 0) {
        step = 0.679 * step;
      }

      direction = direction * (-1);
    }

    COMBIGRID_OUT_LEVEL2( verb , "iterNr=" << iterNr << " , error="
                          << err1 << " , lamdba=" << lamdba1);
    iterNr = iterNr + 1;
  }

  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_FG_crossvalidation    Final initial error="
                        << err1 << " , lamdba=" << lamdba1 << " , step=" << step);
}


void combigrid::RunTikhonov::computeTikhonov_CG_crossvalidation(
  AbstractCombiGrid* combiG,
  double& lambda,
  std::vector<double>& XCoords,
  std::vector<double>& YPoint) {
  int verb = 6;
  int trainPointNr = 3 * static_cast<int>( YPoint.size() ) / 4;
  int testPointNr = static_cast<int>( YPoint.size() ) - trainPointNr;
  int dim = (int)::round((double)XCoords.size() / (double)YPoint.size());

  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_CG_crossvalidation    , trainPointNr="
                        << trainPointNr << " , testPointNr=" << testPointNr);

  std::vector<double> XCoordsTrain(trainPointNr * dim);
  std::vector<double> YPointTrain(trainPointNr);
  std::vector<double> XCoordsTest(testPointNr * dim);
  std::vector<double> YPointTest(testPointNr);

  // copy the test and train values in special arrays
  for (int mp = 0 ; mp < trainPointNr * dim ; mp ++) {
    XCoordsTrain[mp] = XCoords[mp];
  }

  for (int mp = trainPointNr * dim ; mp < (testPointNr + trainPointNr)*dim ; mp ++) {
    XCoordsTest[mp - trainPointNr * dim] = XCoords[mp];
  }

  for (int mp = 0 ; mp < trainPointNr ; mp ++) {
    YPointTrain[mp] = YPoint[mp];
  }

  for (int mp = trainPointNr ; mp < (testPointNr + trainPointNr) ; mp ++) {
    YPointTest[mp - trainPointNr] = YPoint[mp];
  }

  double lamdba1 , lamdba2 = 2e-3 , err1 , err2 , step = (0.9) * lamdba2;
  int direction = 1 , iterNr = 0;
  lamdba1 = lamdba2;
  // we compute the initial error
  computeTikhonov_CG( combiG , lamdba1 , XCoordsTrain, YPointTrain);
  err1 = measureErrorCG( combiG, dim  , XCoordsTest, YPointTest );
  err2 = err1;
  lamdba2 = lamdba1;
  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_CG_crossvalidation    , initial error="
                        << err1 << " , lamdba=" << lamdba1);

  // iterate till we find the solution , this is a simple minimum search algorithm
  // we either look in one direction or into the other, and we change the "step" size so we should hit the minimum
  while ( (iterNr < 15) && (step > 1e-8) ) {
    //
    lamdba1 = lamdba1 - ((double)direction) * step;
    computeTikhonov_CG( combiG , lamdba1 , XCoordsTrain, YPointTrain);
    err1 = measureErrorCG( combiG, dim  , XCoordsTest, YPointTest );
    COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_CG_crossvalidation   try step lambda=" << lamdba1 <<
                          "  error=" << err1 << " , step=" << step);

    // the step was successful
    if (err1 < err2) {
      err2 = err1;
      lamdba2 = lamdba1;

      if ( direction == 1) {
        step = step / 10.0;
      }
    } else {
      // the step leads to a higher value, change direction
      err2 = err1;
      lamdba2 = lamdba1;

      if (direction == 1) {
        step = 3.0 * step;
      }

      if (direction < 0) {
        step = 0.679 * step;
      }

      direction = direction * (-1);
    }

    COMBIGRID_OUT_LEVEL2( verb , "iterNr=" << iterNr << " , error="
                          << err1 << " , lamdba=" << lamdba1);
    iterNr = iterNr + 1;
  }

  COMBIGRID_OUT_LEVEL2( verb , "computeTikhonov_CG_crossvalidation    Final initial error="
                        << err1 << " , lamdba=" << lamdba1);
}


double combigrid::RunTikhonov::measureErrorFG(FullGridD* fg,
    int dim ,
    std::vector<double>& XCoords,
    std::vector<double>& YPoint) {

  double err = 0.0 , fg_val;
  int verb = 6;
  std::vector<double> coord( dim , 0.0 );
  COMBIGRID_OUT_LEVEL2( verb , "measureErrorFG  , dim=" << dim << " , YPoint.size()=" << YPoint.size() << " , XCoords.size()=" << XCoords.size() );

  // for each point, we measure point wise difference
  for (int p = 0 ; p < (int)YPoint.size() ; p++ ) {
    for (int d = 0 ; d < dim ; d++ ) {
      coord[d] = XCoords[p * dim + d];
    }

    fg_val = fg->eval( coord );
    err = err + (YPoint[p] - fg_val) * (YPoint[p] - fg_val);
  }

  err = ::sqrt( err / (double)YPoint.size());
  return err;
}

double combigrid::RunTikhonov::measureErrorCG(AbstractCombiGrid* combiG,
    int dim ,
    std::vector<double>& XCoords,
    std::vector<double>& YPoint) {

  double err = 0.0 , fg_val;
  std::vector<double> coord( dim , 0.0);

  // for each point, we measure point wise difference
  for (int p = 0 ; p < (int)YPoint.size() ; p++ ) {
    for (int d = 0 ; d < dim ; d++ ) {
      coord[d] = XCoords[p * dim + d];
    }

    fg_val = combiG->eval( coord );
    err = err + (YPoint[p] - fg_val) * (YPoint[p] - fg_val);
  }

  err = ::sqrt( err / (double)YPoint.size());
  return err;
}