/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONTESTMODBSPLINE_HPP
#define OPERATIONTESTMODBSPLINE_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

#include "basis/modbspline/modified_bspline_base.hpp"


namespace sg
{
namespace datadriven
{

  /**
   * This class implements OperationTest for a grids with modified bspline basis functions with a certain degree
   *
   * @version $HEAD$
   */
  class OperationTestModBspline : public OperationTest
  {
  public:
    /**
     * Constructor
     *
     * @param storage the grid's base::GridStorage object
     * @param degree the bspline's degree
     */
    OperationTestModBspline(base::GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

    /**
     * Destructor
     */
    virtual ~OperationTestModBspline() {}

    virtual double test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes);
    virtual double testMSE(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& refValues);
    virtual double testWithCharacteristicNumber(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes, base::DataVector& charaNumbers);

  protected:
    /// Pointer to base::GridStorage object
    base::GridStorage* storage;
    /// Mod Bspline Basis object
    base::SModBsplineBase base;
  };

}
}

#endif /* OPERATIONTESTMODBSPLINE_HPP */
