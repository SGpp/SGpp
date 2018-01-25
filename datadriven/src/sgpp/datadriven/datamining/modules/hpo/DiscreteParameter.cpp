/*
 * DiscreteParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#include "DiscreteParameter.hpp"


namespace sgpp {
namespace datadriven {



std::list<ConfigurationBit> DiscreteParameter::makeConfigBits(){
	int i = 1;
	int c = 2;
	while(c <= max-min+1){
		i++;
		c=c*2;
	}
	return HyperParameter::makeConfigBits(i);
}

int DiscreteParameter::getValue(int* configID){
	if(2**bits.size() < max-min+1){
		double v = 0;
		double m = 1;
	  for(auto bit : bits){
	    v = v + m* bit.evaluate(configID);
	    m = m * 2;
	  }
	  return lround(min+((max-min)*(1+v/(m-1.0))/2));
	}else if(2**bits.size() > max-min+1){
		int c = 2**(bits.size()-1);
		int k = 0;
		while(k+c/2+max-min+1 <= 2**bits.size()){
			c = c/2;
			k = k + c;
		}
		// k = k - c;
		int v = min;
		int m = 1;
		// make sure this doesn't break it
		ConfigurationBit last = bits.back();
		bits.pop_back();
		for(auto bit : bits){
			v = v + m*(bit.evaluate(configID)+1)/2;
			m = m * 2;
		}
		bits.push_back(last);
		int nv = v + c*(last.evaluate(configID)+1)/2;
		if(nv >= min+2**(bits.size()-1) && nv <= max){
			v = nv;
		}
		return v;
	}
}


} /* namespace datadriven */
} /* namespace sgpp */
