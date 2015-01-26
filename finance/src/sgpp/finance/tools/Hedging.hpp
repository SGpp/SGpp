/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEDGING_HPP
#define HEDGING_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace finance {

    /**
     * This class implements the calculations of delta and gamma
     * for hedging. They are written into a file including the corresponding
     * option price.
     *
     * For calculating delta and gamma finite difference with sparse
     * grid evaluations are used.
     */
    class Hedging {
      private:
        /// resoluation in hedging area
        size_t m_res;
        /// epsilon used for calculating finite differences
        double m_eps;
        /// Points at which delta and gamma should be calculated, in Cartesian coordinates
        SGPP::base::DataMatrix* m_hedge_points;
        /// is hedging used with log-transformed grids
        bool m_is_log_transformed;

      public:
        /**
         * Constructor
         *
         * @param hedge_area BoundingBox that describes the full-grid area for which the delta and gamma should be calculated. They must be in Cartesian coordinates!
         * @param resolution number of grid points in every dimension
         * @param eps epsilon used for calculating finite differences
         * @param is_log_transformed set to true if hedging is used with log-transformed grids
         */
        Hedging(SGPP::base::BoundingBox& hedge_area, size_t resolution, double eps, bool is_log_transformed);

        /**
         * Destructor
         */
        ~Hedging();

        /**
         * this routine does the actual calculation of delta and gamma based
         * on a sparse grid and its coefficients.
         *
         * @param sparse_grid the sparse grid
         * @param alpha the sparse grid's coefficients
         * @param file_extension some file extension (e.g. numbering) in order to distinguish different outputs that are written
         */
        void calc_hedging(SGPP::base::Grid& sparse_grid, SGPP::base::DataVector alpha, std::string file_extension);
    };

  }

}

#endif /* HEDGING_HPP */
