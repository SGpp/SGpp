/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEOrthoAdapt.hpp
 */

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOnlineDEOrthoAdapt : public DBMatOnlineDE {
 public:
  /**
   * Constructor
   * Builds DBMatOnlineDEOrthoAdapt object from given offline object
   *
   * @param offline The offline object
   * @param beta The initial weighting factor
   */
  explicit DBMatOnlineDEOrthoAdapt(sgpp::datadriven::DBMatOffline& offline, double beta = 0.);

  /**
   * Rank-one updates/downdates the system matrix, based on the sherman-morrison-formula
   *
   * @param adaptPts Points to refine/coarsen
   * @param doCoarsen triggers between refine/coarsen
   */
  void adapt(sgpp::base::DataMatrix& refinePts, std::vector<size_t>& coarsePts);

 protected:
  // matrix, which holds information about refined/coarsened points
  sgpp::base::DataMatrix b_adapt_matrix_;

  // bool, which tells if refinement/coarsening has been performed yet
  bool b_is_refined;
  /**
   * Solves the system (R + lambda*I) * alpha = b, and obtains alpha
   * The solving is done after offline and online phase and looks as follows:
   * (R + lambda*I)^{-1} * b = ((Q*T^{-1}*Q^{t} + B) * b = alpha,
   * where B yields the refine/coarsen information
   *
   * @param b     The right hand side of the system
   * @param do_cv Specifies, if cross-validation should be done
   */
  void solveSLE(sgpp::base::DataVector& b, bool do_cv) override;

 private:
  /**
   * Rank-one updates/downdates the system matrix, based on the sherman-morrison-formula
   * In the current version of the parent functions, the refinePts already are adapted
   * to the configuration parameter lambda.
   * @param refinePts Points to refine
   * @param coarsePts Points to coarsen
   * @param refine    Decides: true for refine, false for coarsen
   */
  void sherman_morrison_adapt(sgpp::base::DataMatrix& refinePts, std::vector<size_t>& coarsePts,
                              bool refine);
};
}  // namespace datadriven
}  // namespace sgpp
