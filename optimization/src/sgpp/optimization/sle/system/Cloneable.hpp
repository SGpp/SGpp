// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SYSTEM_CLONEABLE_HPP
#define SGPP_OPTIMIZATION_SLE_SYSTEM_CLONEABLE_HPP

#include <memory>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/system/System.hpp>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace system {

        /**
         * Abstract class for "cloneable" linear systems.
         * This class is needed in the case that matrix entry lookups are not possible concurrently
         * (e.g. for hierarchisation systems with Clenshaw-Curtis grids).
         */
        class Cloneable : public System {
          public:
            /**
             * Constructor.
             */
            Cloneable() : System() {
            }

            /**
             * Virtual destructor.
             */
            virtual ~Cloneable() {
            }

            /**
             * Pure virtual method for cloning the linear system.
             * It should return a pointer to the cloned object and it's used for parallel computations
             * (e.g. the getMatrixEntry() method might not be thread-safe).
             *
             * @return pointer to cloned object
             */
            virtual Cloneable* clone() = 0;

            /**
             * @return whether this system derives from Cloneable or not (true)
             */
            bool isCloneable() const {
              return true;
            }
        };

      }
    }
  }
}

#endif
