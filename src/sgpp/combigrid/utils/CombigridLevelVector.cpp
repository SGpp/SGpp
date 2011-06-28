/*
 * CombigridLevelVector.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: kowitz_local
 */

#include "CombigridLevelVector.hpp"


namespace combigrid {

CombigridLevelVector::CombigridLevelVector(int dim){
	levelVec_.push_back(std::vector<int>(dim,LEVELMAX));
	coef_.push_back(1);
}

CombigridLevelVector::CombigridLevelVector(std::vector<int> level) {
	levelVec_.push_back(level);
	coef_.push_back(1);
}

CombigridLevelVector::CombigridLevelVector(std::vector<int> level,double coef) {
	levelVec_.push_back(level);
	coef_.push_back(coef);
}
CombigridLevelVector::CombigridLevelVector(std::vector<std::vector<int> > in,std::vector<double> coef){
	levelVec_=in;
	coef_=coef;
	doAddition();
}

CombigridLevelVector& CombigridLevelVector::operator =(const CombigridLevelVector & rhs){
	if (this==&rhs) return (*this);
	levelVec_=rhs.getLevelVec();
	coef_=rhs.getCoef();
	return (*this);
}

const CombigridLevelVector CombigridLevelVector::operator *(const CombigridLevelVector &b) const {
	std::vector< std::vector<int> > result;
	std::vector<int> buffer(b.getDim());
	std::vector<double> c;
	for (int i = 0; i < getN(); ++i) {
		for (int j = 0; j < b.getN(); ++j) {
			for (int k = 0; k < b.getDim(); ++k) {
				buffer[k]=levelVec_[i][k]<b.levelVec_[j][k]? levelVec_[i][k]:b.levelVec_[j][k];
			}
			result.push_back(buffer);
			c.push_back(coef_[i]*b.getCoef()[j]);
		}
	}

	return CombigridLevelVector(result,c);
}
const CombigridLevelVector CombigridLevelVector::operator -(const CombigridLevelVector &b) const {
	CombigridLevelVector result(*this),inVec(b);
	for (int i = 0; i < b.getN(); ++i) {
		inVec.coef_[i]*=-1.0;
	}
	return result+inVec;
}

const CombigridLevelVector CombigridLevelVector::operator +(const CombigridLevelVector &b) const {
	CombigridLevelVector result(*this);
	for (int i = 0; i < b.getN(); ++i) {
		result.levelVec_.insert(result.levelVec_.end(),b.getLevelVec()[i]);
		result.coef_.insert(result.coef_.end(),b.getCoef()[i]);
	}
	result.doAddition();
	return result;
}

void CombigridLevelVector::doAddition(){
	for (int i = 0; i < (int)levelVec_.size(); ++i) {
		for (int j = i+1; j < (int)levelVec_.size(); ++j) {
			bool same=true;
			for (int k = 0; k < getDim(); ++k) {
				if(levelVec_[i][k]!=levelVec_[j][k]){
					same=false;
					break;
				}
			}
			if(same){
				coef_[i]+=coef_[j];
				coef_[j]=0.0;
			}
		}
	}
	for (int i = getN()-1; i > -1; --i) {
		if(coef_[i]==0.0) {
			levelVec_.erase(levelVec_.begin()+i);
			coef_.erase(coef_.begin()+i);
		}
	}
}

void CombigridLevelVector::printLevelVec(){
	for (int i = 0; i < getN(); ++i) {
		for (int j = 0; j < getDim(); ++j) {
			std::cout<<levelVec_[i][j]<<"\t";
		}
		std::cout<<" | "<<coef_[i]<<std::endl;
	}
}

CombigridLevelVector CombigridLevelVector::getCombiLevels(std::vector<CombigridLevelVector> in){
	CombigridLevelVector unity(in[0].getDim());
	CombigridLevelVector erg=unity - in[0];
	for (int i = 1; i < (int)in.size(); ++i) {
		erg=erg*(unity-in[i]);
	}
	erg=unity-erg;
	return erg;
}

CombigridLevelVector CombigridLevelVector::getCombiLevels(std::vector<std::vector<int> > in){
	std::vector<CombigridLevelVector> buffer;
	for (int i = 0; i < (int)in.size(); ++i) {
		buffer.push_back(CombigridLevelVector(in[i]));
	}
	return getCombiLevels(buffer);

}

std::vector<CombigridLevelVector> CombigridLevelVector::split(){
	std::vector<CombigridLevelVector> buffer(0);
	for (int i = 0; i < getN(); ++i) {
		buffer.push_back(CombigridLevelVector(levelVec_[i],coef_[i]));
	}
	return buffer;
}

CombigridLevelVector CombigridLevelVector::getCombiLevels(CombigridLevelVector in){
	in.doAddition();
	for (int i = 0; i < in.getN(); ++i) {
		if(in.getCoef()[i]!=1) return NULL;
	}
	return getCombiLevels(in.split());

}

CombigridLevelVector CombigridLevelVector::getChanges(std::vector<int> in){
	CombigridLevelVector inVector(in);
	CombigridLevelVector unity(in.size());
//	CombigridLevelVector current(levelVec_,coef_);
	std::vector<CombigridLevelVector> currentVec=(*this).split();
	for(unsigned int i=0;i<levelVec_.size();i++){
		inVector=inVector*(unity-currentVec[i]);
	}
	inVector.doAddition();
	return inVector;

}

void CombigridLevelVector::update(std::vector<int> in){
	(*this)=(*this)+(*this).getChanges(in);
}

}


