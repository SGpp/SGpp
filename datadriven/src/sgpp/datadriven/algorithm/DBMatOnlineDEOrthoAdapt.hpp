// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL
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
   * Performs refinement/coarsening on the offline and/or online objects
   * Returns points, which could not be coarsened
   *
   * @param deletedPoints list of indices of last coarsened points
   * @param newPoints amount of refined points
   * @param lambda the regularization parameter
   * @return vector of indices whose points couldn't be coarsened
   */
  std::vector<size_t> adapt(size_t newPoints, std::list<size_t> deletedPoints, double lambda);

  /**
   * Returns the additive component of the sherman-morrison-formula, which
   * yields all the information about the refined points
   */
  sgpp::base::DataMatrix& getB() { return this->b_adapt_matrix_; };

  /**
   * Adds new DataVector to list of refined points
   * For testing purposes only
   *
   * @param the DataVector to add
   */
  void add_new_refine_point(sgpp::base::DataVector& x) { this->refined_points_.push_back(x); };

  std::vector<sgpp::base::DataVector>* getRefinedPointsPointer() {
    return &(this->refined_points_);
  };

  /**
   * Rank-one updates/downdates the system matrix, based on the sherman-morrison-formula
   * In the current version of the function, the refinePts already are adapted
   * to the regularization parameter lambda.
   *
   * @param newPoints number of refined points
   * @param refine decides: true for refine, false for coarsen
   * @param coarsen_indices the indices of points to coarsen sorted in descendant order
   */
  void sherman_morrison_adapt(size_t newPoints, bool refine,
                              std::vector<size_t> coarsen_indices = {});

  // holds all prior refined points, to know what to coarsen later on
  std::vector<sgpp::base::DataVector> refined_points_;


 protected:
  // matrix, which holds information about refined/coarsened points
  sgpp::base::DataMatrix b_adapt_matrix_;


  // points to end of refined_points_, !which already were processed!
  size_t current_refine_index;

  // tells if refinement/coarsening has been performed on the matrix b_adapt_matrix_ yet
  bool b_is_refined;

  /**
   * Solves the system (R + lambda*I) * alpha = b, and obtains alpha
   * The solving is done after offline and online phase and works as follows:
   * (R + lambda*I)^{-1} * b = ((Q*T^{-1}*Q^{t} + B) * b = alpha,
   * where B yields the refine/coarsen information
   *
   * @param b The right hand side of the system
   * @param do_cv Specifies, if cross-validation should be done (todo: currently not implemented)
   */
  void solveSLE(sgpp::base::DataVector& b, bool do_cv) override;

 private:
  /**
   * Computes the L_2 products of the refined gridpoints and pushes them into the
   * refined_points_ container member. The computed vectors of the products correspond
   * to rows/columns of the lhs matrix
   *
   * @param newPoints The number of points to refine
   * @param newLambda The regularization coefficient added to the diagonal elements
   */
  void compute_L2_gridvectors(size_t newPoints, double newLambda);
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* USE_GSL */
