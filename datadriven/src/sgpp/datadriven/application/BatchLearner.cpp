// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/BatchLearner.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <cmath>
#include <deque>
#include <algorithm>
#include <utility>
#include <limits>

// using base;
using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::numeric_limits;
using std::max;
using std::deque;

namespace sgpp {
namespace datadriven {

BatchLearner::BatchLearner(base::BatchConfiguration batchConfig,
                           base::RegularGridConfiguration gridConfig,
                           solver::SLESolverConfiguration solverConfig,
                           base::AdpativityConfiguration adaptConfig) {
  batchConf = batchConfig;
  gridConf = gridConfig;
  solverConf = solverConfig;
  adaptConf = adaptConfig;

  // cofigure solver
  if (solverConf.type_ == solver::SLESolverType::CG)
    myCG = new solver::ConjugateGradients(solverConf.maxIterations_, solverConf.eps_);
  else if (solverConf.type_ == solver::SLESolverType::BiCGSTAB)
    myCG = new solver::BiCGStab(solverConf.maxIterations_, solverConf.eps_);
  else  // not supported
    throw base::application_exception("BatchLearner: An unsupported SLE solver type was chosen!");

  // open file
  reader.open(batchConf.filename_.c_str());

  if (!reader) {
    cout << "ERROR: file does not exist: " << batchConf.filename_ << endl;
    throw 20;
  }
}

// read input, save result to dataFound and last entry to classFound,
// based on ARFFTools.readARFF(..)
void BatchLearner::stringToDataVector(string input, base::DataVector& dataFound, int& classFound) {
  size_t cur_pos = 0;
  size_t cur_find = 0;
  string cur_value;
  double dbl_cur_value;

  base::DataVector temprow(dimensions);

  for (size_t j = 0; j <= dimensions; j++) {
    cur_find = input.find(",", cur_pos);
    cur_value = input.substr(cur_pos, cur_find - cur_pos);
    dbl_cur_value = atof(cur_value.c_str());

    if (j == dimensions)
      classFound = static_cast<int>(dbl_cur_value);
    else
      temprow.set(j, dbl_cur_value);

    cur_pos = cur_find + 1;
  }

  dataFound.resize(temprow.getSize());
  dataFound.copyFrom(temprow);
}

// if (mapdata): data is mapped, dataFound/classesFound will not be changed
void BatchLearner::stringToDataMatrix(string& input, base::DataMatrix& dataFound,
                                      base::DataVector& classesFound, bool mapData) {
  if (mapData) {
    dataInBatch.clear();
  } else {
    dataFound.resize(0, dimensions);
    dataFound.setAll(0.0f);  // needed?
    classesFound.resize(0);
    classesFound.setAll(0.0f);  // needed?
  }

  // as data in input is provided as a string containing all lines: split them at '\n'
  // to vector data to iterate over
  vector<string> data;
  size_t start = 0, end = 0;

  while ((end = input.find('\n', start)) != string::npos) {
    data.push_back(input.substr(start, end - start));
    start = end + 1;
  }

  data.push_back(input.substr(start));

  // iterate over all found lines
  for (size_t i = 0; i < data.size() - 1; i++) {
    size_t classesHere = static_cast<size_t>(count(data[i].begin(), data[i].end(), ','));

    // the first data entry ever found determines the number of dimensions of the data
    if (dimensions == 0) {
      dimensions = classesHere;
      dataFound.resize(dimensions, 0);
      dataFound.setAll(0.0f);
    } else if (classesHere != dimensions) {
      cout << "skipping line " << i
           << " because it contains too few or too many classes (previous/now): " << dimensions
           << "/" << classesHere << endl;
      continue;
    }

    base::DataVector lineData(0);
    int lineClass = -1;
    stringToDataVector(data[i], lineData, lineClass);

    if (mapData) {
      if (dataInBatch.find(lineClass) ==
          dataInBatch.end())  // first data entry for this class in this batch
        dataInBatch.insert(std::pair<int, base::DataMatrix*>(
            lineClass, new base::DataMatrix(0, dimensions, -1.0)));

      // add found data entry to correct base::DataMatrix in map
      dataInBatch.at(lineClass)->appendRow(lineData);
    } else {
      dataFound.appendRow(lineData);
      classesFound.append(static_cast<double>(lineClass));
    }
  }
}

base::DataVector BatchLearner::applyWeight(base::DataVector alpha, int grid) {
  // wMode 4: only the last batch is relevant -> no further calculations needed, alpha = new
  if (batchConf.wMode_ == 4) return alpha;

  if (alphaStorage.find(grid) == alphaStorage.end()) {
    // if there is no previous alpha known for this grid: store this alpha and return the input
    deque<base::DataVector> temp;
    alphaStorage.insert(std::make_pair(grid, temp));
    alphaStorage.at(grid).push_back(alpha);
    return alpha;
  }

  // wMode 5: weigh old alpha with new alpha by occurences
  if (batchConf.wMode_ == 5) {
    double k = static_cast<double>(dataInBatch.at(grid)->getNrows());
    double n = static_cast<double>(occurences.at(grid));
    double wNew = max(k / (n + k), static_cast<double>(batchConf.wArgument_));

    if (batchConf.verbose_) cout << "old weight: " << 1.0 - wNew << " new weight: " << wNew << endl;

    base::DataVector* old = new base::DataVector(alphaStorage.at(grid)[0]);  // old alpha
    old->mult(1.0 - wNew);
    base::DataVector* now = new base::DataVector(alpha);  // new alpha
    now->mult(wNew);
    old->resizeZero(now->getSize());  // refinement
    now->add(*old);
    alphaStorage.at(grid)[0].resizeZero(now->getSize());  // refinement
    alphaStorage.at(grid)[0].copyFrom(*now);
    return *now;
  } else {  // store all alphas if not in wMode 4 or 5
    alphaStorage.at(grid).push_back(alpha);
  }

  // move old alphas to add new one
  if (alphaStorage.at(grid).size() > batchConf.stack_ && batchConf.stack_ != 0) {
    // remove the last item (assumtion: only need to delete one bcs only added here
    alphaStorage.at(grid).pop_front();
  }

  size_t count = alphaStorage.at(grid).size();  // count of old alphas available for calculation
  vector<double> factors;

  // previous alphas exist
  // calc factors
  double sum = 0.0f;

  for (size_t i = 0; i < count; i++) {
    if (batchConf.wMode_ == 0) {
      factors.push_back(1.0);  // temp: all alphas are equal
    } else if (batchConf.wMode_ == 1) {
      factors.push_back(static_cast<double>(i + 1) * batchConf.wArgument_);  // linear
    } else if (batchConf.wMode_ == 2) {
      factors.push_back(pow(batchConf.wArgument_, static_cast<double>(i + 1)));  // exp
    } else if (batchConf.wMode_ == 3) {
      factors.push_back(static_cast<double>(batchConf.wArgument_) /
                        static_cast<double>(i + 1));              // 1/x bzw arg/x
    } else if (batchConf.wMode_ != 4 && batchConf.wMode_ != 5) {  // 4 and 5 treated elsewhere
      cout << "unsupported weighting mode (mode/arg): " << batchConf.wMode_ << "/"
           << batchConf.wArgument_ << endl;
      throw 42;
    }

    sum += factors[i];
  }

  // calc new alphass

  base::DataVector* temp = new base::DataVector(alpha.getSize());
  temp->setAll(0.0f);

  for (size_t i = 0; i < count; i++) {
    alphaStorage.at(grid)[i].resizeZero(
        temp->getSize());  // should have been done already, just to be save
    base::DataVector* t2 = new base::DataVector(alphaStorage.at(grid)[i]);
    // if (batchConf.verbose_)
    //  cout << "mult with " << factors[i] << "/" << sum <<  " = " << factors[i]*1.0f/sum << endl;
    t2->mult(factors[i] * 1.0f / sum);
    temp->add(*t2);
  }

  return *temp;
}

// predict on the data, return the vector of classes predicted
base::DataVector BatchLearner::predict(base::DataMatrix& testDataset, bool updateNorm) {
  if (updateNorm) {
    // update norm factors
    for (auto const& p : grids) {
      //  for each grid
      double evalsum = 0;

      for (size_t x = 0; x < batchConf.samples_; x++) {
        // generate points per grid
        base::DataVector pt(dimensions);

        // only sample within [0.1,0.9]
        for (size_t d = 0; d < dimensions; d++)  // generate point
          pt[d] = 0.1 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (0.9)));

        // add norm factor
        std::unique_ptr<base::OperationEval> opEval(
            op_factory::createOperationEval(*grids.at(p.first)));
        double temp = opEval->eval(*alphaVectors.at(p.first), pt);

        if (batchConf.verbose_ && fabs(temp) > 100)
          cout << "warning abs>100: " << temp << " for " << pt.toString() << endl;

        evalsum += temp;
      }

      evalsum = evalsum / static_cast<double>(batchConf.samples_);
      // update the normFactor
      normFactors.at(p.first) = evalsum;

      if (batchConf.verbose_)
        cout << "grid " << p.first << ", norm factor: " << normFactors.at(p.first) << endl;
    }
  }

  // predict
  base::DataVector result(testDataset.getNrows());
  result.setAll(-2.0f);

  for (unsigned int i = 0; i < testDataset.getNrows(); i++) {
    base::DataVector pt(testDataset.getNcols());
    testDataset.getRow(i, pt);
    // Compute maximum of all density functions:
    int max_index = -1;
    double max = -1.0f * numeric_limits<double>::max();

    for (auto const& g : grids) {
      std::unique_ptr<base::OperationEval> Eval(op_factory::createOperationEval(*g.second));
      // posterior = likelihood*prior
      double res = Eval->eval(*alphaVectors.at(g.first), pt);

      if (batchConf.samples_ != 0) res /= normFactors.at(g.first);

      // this->prior[class_index]
      if (res > max) {
        max = res;
        max_index = g.first;
      }
    }

    result[i] = static_cast<double>(max_index);
  }

  return result;
}

void BatchLearner::processBatch(string workData) {
  base::DataMatrix temp(0, 0);
  base::DataVector temp2(0);
  stringToDataMatrix(workData, temp, temp2, true);

  // data is now mapped to dataInBatch
  // iterate over the found classes
  for (auto const& p : dataInBatch) {
    if (grids.find(p.first) == grids.end()) {
      // this class has never been seen before
      // init new grid etc
      grids.insert(std::pair<int, base::LinearGrid*>(p.first, new base::LinearGrid(dimensions)));
      occurences.insert(std::pair<int, int>(p.first, 0));
      // Generate regular Grid with LEVELS Levels
      grids.at(p.first)->getGenerator().regular(gridConf.level_);

      if (batchConf.verbose_)
        cout << "found new class " << p.first << ", points in grid " << p.first << ": "
             << grids.at(p.first)->getSize() << endl;

      alphaVectors.insert(std::pair<int, base::DataVector*>(
          p.first, new base::DataVector(grids.at(p.first)->getSize())));
      alphaVectors.at(p.first)->setAll(0.0);
      normFactors.insert(std::pair<int, double>(p.first, 1));
    }

    // calculate the new surplusses for the items found of this class in this batch
    // Solve the system for every class and store coefficients in newAlpha
    // set up everything to be able to solve
    base::DataVector newAlpha(grids.at(p.first)->getSize());
    newAlpha.setAll(0.0);
    base::OperationMatrix* id = op_factory::createOperationIdentity(*grids.at(p.first));
    datadriven::DensitySystemMatrix DMatrix(*grids.at(p.first), *dataInBatch.at(p.first), id,
                                            batchConf.lambda_);
    base::DataVector rhs(grids.at(p.first)->getSize());
    DMatrix.generateb(rhs);
    myCG->setMaxIterations(solverConf.maxIterations_);
    myCG->setEpsilon(solverConf.eps_);
    // solve euqation to get new alpha
    myCG->solve(DMatrix, newAlpha, rhs, false, false, -1.0);

    // apply weighting
    alphaVectors.at(p.first)->copyFrom(applyWeight(newAlpha, p.first));

    // refine the grids
    if (batchConf.refineEvery_ != 0 && batchnum != 0 && batchnum % batchConf.refineEvery_ == 0) {
      if (batchConf.verbose_) cout << "refining ..." << endl;

      base::SurplusRefinementFunctor myRefineFunc(*alphaVectors.at(p.first), adaptConf.noPoints_,
                                                  adaptConf.threshold_);
      grids.at(p.first)->getGenerator().refine(myRefineFunc);
      // change alpha, zeroes to new entries until they will be filled
      alphaVectors.at(p.first)->resizeZero(grids.at(p.first)->getSize());
    }

    occurences.at(p.first) += p.second->getSize();
  }
}

void BatchLearner::trainBatch() {
  int startLine = dataLine;
  string line = "";
  string test = "";
  unsigned int ts = 0;

  // collect data from arff to process
  while (bs + ts != batchConf.batchsize_ + batchConf.testsize_) {
    isFinished = reader.eof();

    if (isFinished) {
      reader.close();
      return;
    }

    // read another line from arff
    getline(reader, line);

    if (!reachedData && line.find("@DATA") != string::npos) {
      reachedData = true;
      continue;
    }

    if (!reachedData) continue;

    // only count lines after @DATA
    dataLine++;

    if (bs < batchConf.batchsize_) {
      batch += line + "\n";
      bs++;
    } else if (ts < batchConf.testsize_) {
      test += line + "\n";
      ts++;
    }
  }

  // process the batch
  if (batchConf.verbose_)
    cout << "Processing Data " << batchnum << " @" << startLine << "-" << dataLine << " ..."
         << endl;

  processBatch(batch);

  // automatic testing if wanted by user
  if (batchConf.testsize_ > 0 && test.size() > 0) {
    base::DataMatrix testData(0, 0);
    base::DataVector testClasses(0);
    stringToDataMatrix(test, testData, testClasses, false);
    // predict
    base::DataVector result = predict(testData, batchConf.samples_ != 0);

    if (batchConf.verbose_) cout << "Testing Data (" << result.getSize() << " items)..." << endl;

    // caluclate accuracy
    if (result.getSize() > 0) {
      if (batchConf.verbose_) {
        cout << "result: " << result.toString() << endl;
        cout << "should: " << testClasses.toString() << endl;
      }

      // count correct entries
      int correct = 0;

      for (size_t i = 0; i < result.getSize(); i++) {
        if (result.get(i) == testClasses.get(i)) correct++;
      }

      // calc accuracy for this batch and all tests
      t_total += static_cast<int>(result.getSize());
      t_correct += correct;
      acc_current = 100.0 * correct / static_cast<double>(result.getSize());
      acc_global = 100.0 * t_correct / static_cast<double>(t_total);
      // output accuracy
      cout << "batch:\t" << acc_current << "% (" << correct << "/" << result.getSize() << ")"
           << endl;
      cout << "total:\t" << acc_global << "% (" << t_correct << "/" << t_total << ")" << endl;
    }
  }

  batch = "";
  batch += test;
  bs = ts;
  batchnum++;
}

}  // namespace datadriven
}  // namespace sgpp
