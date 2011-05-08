/*
 * CombigridLevelVector.hpp
 *
 *  Created on: Apr 28, 2011
 *      Author: kowitz_local
 */

#ifndef COMBIGRIDLEVELVECTOR_HPP_
#define COMBIGRIDLEVELVECTOR_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"

namespace combigrid {

class CombigridLevelVector {
public:
	CombigridLevelVector(std::vector<int> level);
	CombigridLevelVector(int dim);
	CombigridLevelVector(){levelVec_.resize(0);coef_.resize(0);};
	CombigridLevelVector(std::vector<int> level, double coef);
	CombigridLevelVector(std::vector<std::vector<int> > in,
			std::vector<double> coef);

	std::vector<std::vector<int> > getLevelVec() const {
		return levelVec_;
	}
	std::vector<int> getLevelVecSingle(int i) const {
		return levelVec_[i];
	}
	std::vector<double> getCoef() const {
		return coef_;
	}

	int getDim() const {
		return levelVec_[0].size();
	}

	int getN() const {
		return levelVec_.size();
	}

	CombigridLevelVector& operator=(const CombigridLevelVector & rhs);
	const CombigridLevelVector operator*(const CombigridLevelVector & b) const;
	const CombigridLevelVector operator+(const CombigridLevelVector & b) const;
	const CombigridLevelVector operator-(const CombigridLevelVector & b) const;

	void doAddition();
	void printLevelVec();

	std::vector<CombigridLevelVector> split();

	/**
	 * Function creating the levels of a combigrid containing all subgrids
	 * provided in the Combigridlevelvectors
	 */
	static CombigridLevelVector getCombiLevels(std::vector<CombigridLevelVector> in);
	static CombigridLevelVector getCombiLevels(std::vector<std::vector<int> > in);
	static CombigridLevelVector getCombiLevels(CombigridLevelVector in);

protected:
	std::vector<std::vector<int> > levelVec_;
	std::vector<double> coef_;
	static const int LEVELMAX = 128;
};

}

#endif /* COMBIGRIDLEVELVECTOR_HPP_ */
