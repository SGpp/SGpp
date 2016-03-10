// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef UPDOWNFOUROPDIMS_HPP
#define UPDOWNFOUROPDIMS_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

#include <sgpp/globaldef.hpp>

#include <vector>
#include <map>

namespace sgpp {
namespace pde {

/**
 * Implements the Up/Down scheme with four dimensions with special operations: i,j,k,l
 *
 */
class UpDownFourOpDims : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   * @param coef 4d tensor that contains the constant coefficients of this operation
   */
  UpDownFourOpDims(sgpp::base::GridStorage* storage, double***** coef);

  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit UpDownFourOpDims(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~UpDownFourOpDims();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 protected:
  typedef sgpp::base::GridStorage::grid_iterator grid_iterator;

  /// Function pointer type. This is used in fnMap to map the particular dimension situation to the
  /// relevant method handler.
  typedef void (sgpp::pde::UpDownFourOpDims::*MFP)(sgpp::base::DataVector&, sgpp::base::DataVector&,
                                                   size_t, size_t, size_t, size_t, size_t);

  /// Pointer to the grid's storage object
  sgpp::base::GridStorage* storage;
  /// Pointer to the coefficients of this bilinear form
  double**** coefs;
  /// algorithmic dimensions, operator is applied in this dimensions
  const std::vector<size_t> algoDims;
  /// number of algorithmic dimensions
  const size_t numAlgoDims_;
  /// max number of parallel stages (dimension recursive calls)
  static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

  /// Map of integer to function pointer. This is used to map the dimension situation to the
  /// relevant method handler.
  std::map<size_t, MFP> fnMap;

  /**
   * Utility method to generate the fnMap member for mappings.
   */
  void generateMap();

  /**
   * Recursive procedure for updown, parallel version using OpenMP 3
   *
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   * @param dim the current dimension
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   * @param op_dim_three the dimension in which to use the third gradient
   * @param op_dim_four the dimension in which to use the fourth gradient
   */
  void updown(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
              size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

  // Unidirectional
  void specialOpUnidirectional(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                               size_t dim, size_t op_dim_one, size_t op_dim_two,
                               size_t op_dim_three, size_t op_dim_four);

  // Singles
  void specialOpOne(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                    size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
  void specialOpTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                    size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
  void specialOpThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                      size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                      size_t op_dim_four);
  void specialOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim,
                     size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
  void specialOpX(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                  void (sgpp::pde::UpDownFourOpDims::*pt2UpFunc)(sgpp::base::DataVector&,
                                                                 sgpp::base::DataVector&, size_t),
                  void (sgpp::pde::UpDownFourOpDims::*pt2DownFunc)(sgpp::base::DataVector&,
                                                                   sgpp::base::DataVector&, size_t),
                  size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                  size_t op_dim_four);

  // Doubles

  /**
   * If the current dimension is equal to the both special operation dimensions one and two.
   *
   * @param alpha the coefficients of the grid points
   * @param result the result of the operations
   * @param dim the current dimension in the recursion
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   * @param op_dim_three the dimension in which to use the third gradient
   * @param op_dim_four the dimension in which to use the fourth gradient
   */
  void specialOpOneAndOpTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                            size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                            size_t op_dim_four);

  /**
   * If the current dimension is equal to the both special operation dimensions one and three.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                              size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                              size_t op_dim_four);

  /**
   * If the current dimension is equal to the both special operation dimensions one and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                             size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                             size_t op_dim_four);

  /**
   * If the current dimension is equal to the both special operation dimensions two and three.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpTwoAndOpThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                              size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                              size_t op_dim_four);

  /**
   * If the current dimension is equal to the both special operation dimensions two and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpTwoAndOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                             size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                             size_t op_dim_four);

  /**
   * If the current dimension is equal to the both special operation dimensions three and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpThreeAndOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                               size_t dim, size_t op_dim_one, size_t op_dim_two,
                               size_t op_dim_three, size_t op_dim_four);

  /**
   * If the current dimension is equal to the all special operation dimensions one, two and three.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpTwoAndOpThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim, size_t op_dim_one, size_t op_dim_two,
                                      size_t op_dim_three, size_t op_dim_four);

  /**
   * If the current dimension is equal to the all special operation dimensions one, two and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpTwoAndOpFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim, size_t op_dim_one, size_t op_dim_two,
                                     size_t op_dim_three, size_t op_dim_four);

  /**
   * If the current dimension is equal to the all special operation dimensions one, three and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim,
                                       size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                                       size_t op_dim_four);

  /**
   * If the current dimension is equal to the all special operation dimensions two, three and four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpTwoAndOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim,
                                       size_t op_dim_one, size_t op_dim_two, size_t op_dim_three,
                                       size_t op_dim_four);

  /**
   * If the current dimension is equal to the all special operation dimensions one, two, three and
   * four.
   * For an explanation of the parameters of this method, see the documentation for the method
   * specialOpOneAndOpTwo in this class.
   */
  void specialOpOneAndOpTwoAndOpThreeAndOpFour(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataVector& result, size_t dim,
                                               size_t op_dim_one, size_t op_dim_two,
                                               size_t op_dim_three, size_t op_dim_four);

  /**
   * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the up-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   * @param dim dimension in which to apply the up-part
   */
  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the down-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   * @param dim dimension in which to apply the down-part
   */
  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i.
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimOne(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                            size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOne(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                          size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to j.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                            size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to j.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                          size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                              size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimThree(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                            size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                             size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                           size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and j.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and j.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                     size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimThree(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to j and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to j and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to j and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to j and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                          sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                        sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and j and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and j and k.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimTwoAndOpDimThree(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and j and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and j and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimTwoAndOpDimFour(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result,
                                                     size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to j and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                     sgpp::base::DataVector& result,
                                                     size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to j and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                   sgpp::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i and j and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                                sgpp::base::DataVector& result,
                                                                size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i and j and k and l.
   * For an explanation of the parameters of this method, see the documentation for the method
   * downOpDimOne in this class.
   */
  virtual void upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sgpp::base::DataVector& alpha,
                                                              sgpp::base::DataVector& result,
                                                              size_t dim) = 0;
};
}  // namespace pde
}  // namespace sgpp

#endif /* UPDOWNFOUROPDIMS_HPP */
