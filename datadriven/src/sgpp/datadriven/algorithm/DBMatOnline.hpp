// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef DBMATONLINE_HPP_
#define DBMATONLINE_HPP_

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

/**
 * Class for objects that can be used in the online step of the classification
 * (The classification is divided into two parts: the offline step that does not
 * depend on the actual data and the online step that depends on the data)
 */
class DBMatOnline {
 public:
  /**
   * Constructor
   */
  DBMatOnline();

  /**
   * Constructor
   *
   * @param o a offline object
   */
  DBMatOnline(DBMatOffline* o);

  /**
   * Destructor
   */
  virtual ~DBMatOnline();

  /**
   * Reads an offline object
   *
   * @param o the offline object
   */
  virtual void readOffline(DBMatOffline* o);

  /**
   * Changes the weighting factor for the regularization term, 
   * if possible (might depend on the kind of decomposition for classification)
   */
  virtual void setLambda(double lambda);

  /**
   * Returns a pointer to the offline object
   */
  DBMatOffline* getOffline();

 protected:
  DBMatOffline* offlineObject_;
};

#endif /* DBMATONLINE_HPP_ */

#endif /* USE_GSL */

