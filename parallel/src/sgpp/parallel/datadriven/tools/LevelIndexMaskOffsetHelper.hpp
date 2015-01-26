/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef LEVELINDEXMASKOFFSETHELPER_HPP
#define LEVELINDEXMASKOFFSETHELPER_HPP

#include <sgpp/parallel/datadriven/basis/common/KernelBase.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    // comment is true for both DP and SP
    // use constructor to do actual work, in code where we use this, it looks like a
    // function call (and we don't have to return anything, so this is nicer)
    // C must have the following members:
    // level_
    // index_
    // mask_
    // offset_
    // storage_
    //
    // if these members are not visible, add the following lines (for each KernelType to
    // support one line) to the class definition of the classes (NameOfClass) you want
    // to use this helper with. When you have all the fields defined in a superclass, using
    // this superclass as template parameter (and only adding the friend struct ... there)
    // avoids writing too much boilerplate code.
    // friend struct LevelIndexMaskOffsetHelper::rebuild<Mask, NameOfClass>;
    // friend struct LevelIndexMaskOffsetHelper::rebuild<Standard, NameOfClass>;

    namespace LevelIndexMaskOffsetHelper {

      template<KernelType T, typename C> struct rebuild {
        rebuild(C* op);
      };
      template<typename C> struct rebuild<Standard, C> {
        rebuild(C* op);
      };
      template<typename C> struct rebuild<Mask, C> {
        rebuild(C* op);
      };

      template<KernelType T, typename C>  rebuild<T, C>::rebuild(C* op) {
        throw base::operation_exception("The rebuild operation for the specified KernelType is not implemented.");
      }
      template<typename C>
      inline rebuild<Standard, C>::rebuild(C* op) {
        if (op->level_ != NULL)
          delete op->level_;

        if (op->index_ != NULL)
          delete op->index_;

        op->level_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());
        op->index_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());

        op->storage_->getLevelIndexArraysForEval(*(op->level_), *(op->index_));
      }
      template<typename C>
      inline rebuild<Mask, C>::rebuild(C* op) {
        if (op->level_ != NULL)
          delete op->level_;

        if (op->index_ != NULL)
          delete op->index_;

        if (op->mask_ != NULL)
          delete op->mask_;

        if (op->offset_ != NULL)
          delete op->offset_;

        op->level_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());
        op->index_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());
        op->mask_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());
        op->offset_ = new SGPP::base::DataMatrix(op->storage_->size(), op->storage_->dim());

        op->storage_->getLevelIndexMaskArraysForModEval(*(op->level_), *(op->index_), *(op->mask_), *(op->offset_));
      }
    }


    namespace LevelIndexMaskOffsetHelperSP {
      template<KernelType T, typename C> struct rebuild {
        rebuild(C* op);
      };
      template<typename C> struct rebuild<Standard, C> {
        rebuild(C* op);
      };
      template<typename C> struct rebuild<Mask, C> {
        rebuild(C* op);
      };

      template<KernelType T, typename C>  rebuild<T, C>::rebuild(C* op) {
        throw base::operation_exception("The rebuild operation for the specified KernelType is not implemented.");
      }
      template<typename C>
      inline rebuild<Standard, C>::rebuild(C* op) {
        if (op->level_ != NULL)
          delete op->level_;

        if (op->index_ != NULL)
          delete op->index_;

        op->level_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());
        op->index_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());

        op->storage_->getLevelIndexArraysForEval(*(op->level_), *(op->index_));
      }
      template<typename C>
      inline rebuild<Mask, C>::rebuild(C* op) {
        if (op->level_ != NULL)
          delete op->level_;

        if (op->index_ != NULL)
          delete op->index_;

        if (op->mask_ != NULL)
          delete op->mask_;

        if (op->offset_ != NULL)
          delete op->offset_;

        op->level_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());
        op->index_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());
        op->mask_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());
        op->offset_ = new SGPP::base::DataMatrixSP(op->storage_->size(), op->storage_->dim());

        op->storage_->getLevelIndexMaskArraysForModEval(*(op->level_), *(op->index_), *(op->mask_), *(op->offset_));
      }
    }

  }
}

#endif // LEVELINDEXMASKOFFSETHELPER_HPP
