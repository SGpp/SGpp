#include "OperationMatrixLTwoDotExplicitLinearBoundary.hpp"
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/data_exception.hpp>

namespace sg {
  namespace pde {

    OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(sg::base::DataMatrix* m,
        sg::base::Grid* grid) :
      ownsMatrix_(false) {
      m_ = m;
      buildMatrix(grid);
    }

    OperationMatrixLTwoDotExplicitLinearBoundary::OperationMatrixLTwoDotExplicitLinearBoundary(sg::base::Grid* grid) :
      ownsMatrix_(true) {
      m_ = new sg::base::DataMatrix(grid->getStorage()->size(),
                                    grid->getStorage()->size());
      buildMatrix(grid);
    }

    OperationMatrixLTwoDotExplicitLinearBoundary::~OperationMatrixLTwoDotExplicitLinearBoundary() {
      if (ownsMatrix_)
        delete m_;
    }

    void OperationMatrixLTwoDotExplicitLinearBoundary::buildMatrix(sg::base::Grid* grid) {
      //Build matrix (in the moment just by multiplying the OperationMatrix with the unit vectors):
      OperationMatrix* opMatrix = sg::op_factory::createOperationLTwoDotProduct(
                                    *grid);

      size_t size = grid->getStorage()->size();
      sg::base::DataVector unit(size);
      unit.setAll(0.0);
      sg::base::DataVector result(size);

      for (size_t i = 0; i < size; i++) {
        //Compute i-th unit vector
        if (i > 0)
          unit.set(i - 1, 0.0);

        unit.set(i, 1.0);

        //Multiply with operation matrix
        opMatrix->mult(unit, result);
        m_->setColumn(i, result);
      }
    }

    void OperationMatrixLTwoDotExplicitLinearBoundary::mult(sg::base::DataVector& alpha,
        sg::base::DataVector& result) {

      size_t nrows = m_->getNrows();
      size_t ncols = m_->getNcols();

      if (alpha.getSize() != ncols || result.getSize() != nrows) {
        throw sg::base::data_exception("Dimensions do not match!");
      }

      double* data = m_->getPointer();

      //Standard matrix multiplication:
      double temp = 0.;
      size_t acc = 0;

      for (size_t i = 0; i < nrows; i++) {
        for (size_t j = 0; j < ncols; j++) {
          temp += data[j + acc] * alpha[j];
        }

        result[i] = temp;
        temp = 0.;
        acc += ncols;
      }
    }

  } /* namespace pde */
} /* namespace sg */
