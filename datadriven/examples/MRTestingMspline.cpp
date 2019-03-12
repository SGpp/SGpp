// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/activeSubspaces/MSplineBasis.hpp>

#include <iostream>

int main() {
  sgpp::base::SGppStopwatch watch;
  sgpp::base::DataVector xi(4);
  for (size_t i = 0; i < xi.getSize(); i++) {
    xi[i] = static_cast<double>(i);
  }
  sgpp::datadriven::MSplineBasis bas(xi);
  size_t degree = 3;
  size_t index = 0;
  double x = 0.5;
  double evalS = 0.0;
  double trevalS = 0.0;
  watch.start();
  for (size_t i = 0; i < 100000000; i++) {
    evalS += bas.eval(degree, index, x);
  }
  double evalTime = watch.stop();
  watch.start();
  for (size_t i = 0; i < 100000000; i++) {
    trevalS += bas.evalTruncated(x);
  }
  double trevalTime = watch.stop();
  std::cout << evalS << "\n";
  std::cout << trevalS << "\n";
  std::cout << evalTime << "\n";
  std::cout << trevalTime << "\n";
}
