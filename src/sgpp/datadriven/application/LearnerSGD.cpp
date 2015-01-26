#include "base/operation/BaseOpFactory.hpp"
#include "LearnerSGD.hpp"
#include <cmath>

 
namespace sg {

  namespace datadriven {

    LearnerSGD::LearnerSGD(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool isVerbose) : Learner(regularization, isRegression, isVerbose) {
    }

	void LearnerSGD::train(
			sg::base::DataMatrix& trainDataset, 
			sg::base::DataVector& classes, 
			sg::base::RegularGridConfiguration& GridConfig, 
		    size_t maxIterations,	
			double eps,
            double lambda,
			double gamma
			) {
		using namespace sg::base;

		// Initialize Grid
		InitializeGrid(GridConfig);
		if (grid_ == NULL)
	      return;

		alpha_->setAll(0.0);

		size_t num_coeff = grid_->getStorage()->size();
		size_t dim = trainDataset.getNcols();

		// Execute SGD
		size_t numIterations = 0;

		while (numIterations < maxIterations)
		{
			// Get random x and y pair
			int k = getRandom((int) trainDataset.getNrows() - 1); 
			DataVector x(dim);
		   	trainDataset.getRow((size_t)k, x);
			double y = classes[k];

			// Calculate delta^n according to [Maier BA, 5.10]:

			// tmp1 = (b_k^T * a^n - y_k) where
			// b_k = (phi_1(x_k) ... phi_N(x_k))
			double tmp1 = grid_->eval(*alpha_, x) - y;

			// delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
			DataVector delta(num_coeff);

			for (unsigned int i=0; i < num_coeff; i++)
			{
				DataVector unit_alpha(num_coeff);
				unit_alpha.setAll(0.0);
				unit_alpha[i] = 1;

				delta[i] = 2 * gamma * ( grid_->eval(unit_alpha, x) * tmp1 + lambda * (*alpha_)[i]  );
			}
			
			// update alpha
			// a^{n+1} = a^n - delta^n 
			for(unsigned int i=0; i < num_coeff; i++)
			{
				(*alpha_)[i] = (*alpha_)[i] - delta[i];
			}

			// check if below eps
			bool is_below_eps = true;
			for (unsigned int i=0; i < num_coeff; i++)
			{
				if (fabs(delta[i]) >= eps)
				{
					is_below_eps = false;
					break;
				}
			}

			if (is_below_eps)
			{
				return;
			}

			numIterations++;
		}

		isTrained_ = true;
	}

	int LearnerSGD::getRandom(int limit)
	{
		int divisor = RAND_MAX/(limit+1);
		int r;

		do { 
			r = rand() / divisor;
		} while (r > limit);

		return r;
	}

	sg::base::DataVector* LearnerSGD::getAlpha()
	{
		return alpha_;
	}

	sg::base::Grid* LearnerSGD::getGrid()
	{
		return grid_;
	}

	LearnerSGD::~LearnerSGD() {
	}
  }
}
