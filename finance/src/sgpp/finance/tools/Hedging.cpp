// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/tools/Hedging.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/EvalCuboidGenerator.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sstream>
#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace finance {

    Hedging::Hedging(SGPP::base::BoundingBox& hedge_area, size_t resolution, float_t eps, bool is_log_transformed) :
      m_res(resolution), m_eps(eps), m_hedge_points(new SGPP::base::DataMatrix(1, hedge_area.getDimensions())), m_is_log_transformed(is_log_transformed) {
      SGPP::base::EvalCuboidGenerator* myEval = new SGPP::base::EvalCuboidGenerator();
      myEval->getEvaluationCuboid(*m_hedge_points, hedge_area, resolution);
      delete myEval;
    }

    Hedging::~Hedging() {
      delete m_hedge_points;
    }

    void Hedging::calc_hedging(SGPP::base::Grid& sparse_grid, SGPP::base::DataVector alpha, std::string file_extension) {
      std::stringstream sfilename;
      sfilename << "hedging_" << file_extension << ".out";

      std::ofstream fileout;
      fileout.open(sfilename.str().c_str());

      SGPP::base::OperationEval* myEval = SGPP::op_factory::createOperationEval(sparse_grid);

      // loop overall evaluation points
      for (size_t i = 0; i < m_hedge_points->getNrows(); i++) {
        base::DataVector curPoint(m_hedge_points->getNcols());
        base::DataVector left(m_hedge_points->getNcols());
        base::DataVector right(m_hedge_points->getNcols());

        m_hedge_points->getRow(i, curPoint);


        // print coordinates into file
        for (size_t j = 0; j < m_hedge_points->getNcols(); j++) {
          fileout << curPoint.get(j) << " ";
        }

        if (m_is_log_transformed == true) {
          for (size_t j = 0; j < m_hedge_points->getNcols(); j++) {
            curPoint.set(j, log(curPoint.get(j)));
          }
        }

        // get option value
        float_t value = myEval->eval(alpha, curPoint);
        fileout << value << " ";

        // calculate derivatives
        float_t delta = 0.0;
        float_t gamma = 0.0;

        for (size_t j = 0; j < m_hedge_points->getNcols(); j++) {
          m_hedge_points->getRow(i, left);
          m_hedge_points->getRow(i, right);

          if (m_is_log_transformed == true) {
            left.set(j, log((left.get(j) - m_eps)));
            right.set(j, log((right.get(j) + m_eps)));
          } else {
            left.set(j, (left.get(j) - m_eps));
            right.set(j, (right.get(j) + m_eps));
          }

          float_t tmp_delta = ((myEval->eval(alpha, right) - myEval->eval(alpha, left)) / (2.0 * m_eps));
          float_t tmp_gamma = ((myEval->eval(alpha, right) - (2.0 * value) + myEval->eval(alpha, left)) / (m_eps * m_eps));

          fileout << tmp_delta << " " << tmp_gamma << " ";

          delta += tmp_delta;
          gamma += tmp_gamma;
        }

        fileout << delta << " " << gamma << std::endl;

        if ((i + 1) % m_res == 0) {
          fileout << std::endl;
        }
      }

      delete myEval;
      fileout.close();
    }

  }

}