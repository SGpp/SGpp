#ifndef INTEGRATION_HPP_
#define INTEGRATION_HPP_

#include <queue>
#include <vector>
#include <functional>
#include <algorithm>

#include <iostream> 

#include <math.h>

template<typename IT, typename DT>
struct RefineNode
{
public:
	RefineNode(DT alpha, const IT* index) : alpha(alpha), index(index)
	{
	}
	
    bool operator<( const RefineNode &n ) const {
        return alpha < n.alpha;
    }

    bool operator>( const RefineNode &n ) const {
        return alpha > n.alpha;
    }
	
	
	DT alpha;
	const IT* index;	
};

template<typename STORAGE, typename DATA, typename FUNC>
class LASGridIntegrator
{
public:
	typedef typename STORAGE::index_type index_type;
	typedef typename STORAGE::index_pointer index_pointer;
	typedef RefineNode<index_type, double> refine_type;
	
	LASGridIntegrator(STORAGE &storage, FUNC functor, double move) : storage(storage), functor(functor), interm(), alpha(), refineQueue(), integral_value(0.0), move(move)
	{
		interm.resize(storage.dim(), DATA());
		
		index_type index;
		for(int i = 0; i < storage.dim(); i++)
		{
			index.push(i, 1, 1);
		}		
		index.rehash();
		
		
		index_type* insert = storage.create(index);
		storage.store(insert);
		
		double value = functor(index);
		for(int i = 0; i < storage.dim(); i++)
		{
			interm[i].push_back(value);
		}
		
		alpha.push_back(value);
		
		refineQueue.push( refine_type( fabs(value), insert ));
		integral_value += value * insert->volume();
		
		
		//createGridPoint(index);
	}
	
	void refine()
	{
		while(refineQueue.size() > 0)
		{ 
			refine_type rNode = refineQueue.top();
			refineQueue.pop();
			
			index_type index(rNode.index);
			
			bool refined = false;
			 			
			for(int d = 0; d < storage.dim(); d++)
			{
				typename index_type::level_type srcLevel;
				typename index_type::index_type srcIndex;
	
				index.get(d, srcLevel, srcIndex);
				
				index.set(d, srcLevel+1, srcIndex*2-1);
				if(!storage.has_key(&index))
				{
					createGridPoint(index);
					refined = true;	
				}
				
				index.set(d, srcLevel+1, srcIndex*2+1);
				if(!storage.has_key(&index))
				{
					createGridPoint(index);
					refined = true;	
				}
				
				index.set(d, srcLevel, srcIndex);
			}
			if(refined)
			{
				break;
			}
		}				
	}
	
	void toStringAlpha(std::ostream &stream)
	{
		for(unsigned int i = 0; i < alpha.size(); i++)
		{
			stream << alpha[i] << ", ";
		}
		stream << std::endl;
	}
	
	void printGrid(std::ostream &stream)
	{
		for(unsigned int i = 0; i < alpha.size(); i++)
		{
			storage[i]->printGrid(stream);
			stream << alpha[i] << std::endl;
		}
		
	}
	
	double volume()
	{
		/*
		typename STORAGE::grid_iterator iter = storage.begin();
		typename STORAGE::grid_const_iterator end = storage.end();
		
		double vol = 0.0;
		
		while(end != iter)
		{
			vol += iter->first->volume() * alpha[iter->second];
			iter++;
		}
		return vol;
		*/
		return integral_value;
	}

protected:
	
	unsigned int createGridPoint(index_type &index)
	{
		index_type* insert = storage.create(index);
		unsigned int seq = storage.store(insert);
		
		
		for(int d = 0; d < storage.dim(); d++)
		{
			interm[d].push_back(0.0);
		}

		double value = 0.0;
		
		for(int d = 0; d < storage.dim(); d++)
		{
			double interpolated = 0.0;
						
			typename index_type::level_type srcLevel;
			typename index_type::index_type srcIndex;
	
			index_pointer left_index = NULL;
			index_pointer right_index = NULL;
			
			double left = 0.0;
			double right = 0.0;
		
			interm[d][seq] = value;
			
			index.get(d, srcLevel, srcIndex);
			{
				typename index_type::level_type wLevel = srcLevel;
				typename index_type::index_type wIndex = srcIndex - 1; 
	
	            while(!(wIndex & 0x01) && wLevel-- > 1)
	            {
	                wIndex = wIndex >> 1;
	            }
	            if (wLevel > 0)
	            {
	                index.set(d, wLevel, wIndex);

					typename STORAGE::grid_iterator iter = storage.find(&index);
					unsigned int seqParent;
					if (storage.end() == iter)
					{
						seqParent = createGridPoint(index);
					}
					else
					{
						seqParent = iter->second;
					}
					left = interm[d][seqParent];
					left_index = storage[seqParent];
					//interpolated += interm[d][seqParent] * storage[seqParent]->influence(d, index);
					//interpolated += interm[d][seqParent] * 0.5;
	            }
			}

			{
				typename index_type::level_type wLevel = srcLevel;
				typename index_type::index_type wIndex = srcIndex + 1; 
	
	            while(!(wIndex & 0x01) && wLevel-- > 1)
	            {
	                wIndex = wIndex >> 1;
	            }
	            if (wLevel > 0)
	            {
	                index.set(d, wLevel, wIndex);

					typename STORAGE::grid_iterator iter = storage.find(&index); 
					unsigned int seqParent;
					if (storage.end() == iter)
					{
						seqParent = createGridPoint(index);
					}
					else
					{
						seqParent = iter->second;
					}
					right = interm[d][seqParent];
					right_index = storage[seqParent];
					//interpolated += interm[d][seqParent] * storage[seqParent]->influence(d, index);
					//interpolated += interm[d][seqParent] * 0.5;
	            }
			}
			
			//totally unportable stuff
			
			double base = 0.0;
			double p = 0.0;
			if(right_index != NULL)
			{
				base = right_index->getAbs(d);
				p += right_index->getAbs(d);
			} 
			else 
			{
				base = 1.0;
				p += 1.0;
			}

			if(left_index != NULL)
			{
				base -= left_index->getAbs(d);
				p += left_index->getAbs(d);
			}
			
			p *= 0.5;
			if(p > 0.5)
			{
				p -= base*(0.5-move);
				interpolated += (1.0-move) * left + move * right;
			}
			else if(p < 0.5)
			{
				p += base*(0.5-move);
				interpolated += move * left + (1.0-move) * right;
			}
			else
			{
				interpolated += 0.5 * (left + right);
			}
			
			insert->move(d, p, base);			
			
			// *******
			//interpolated += (left + right) * 0.5;
			
			index.set(d, srcLevel, srcIndex);
			value -= interpolated;
		}
		
		double func = functor(*insert);
		value += func;
		for(int d = 0; d < storage.dim(); d++)
		{
			interm[d][seq] += func;
		}
		
		//std::cout << value << std::endl;
		alpha.push_back(value);
		refineQueue.push( refine_type( fabs(value), insert ));
		integral_value += value * insert->volume();
		
		return seq;
	}
	
	
	
	
private:
	STORAGE& storage;
	FUNC functor;
	std::vector<DATA> interm;
	DATA alpha;
	std::priority_queue<refine_type, std::vector<refine_type>, std::less<refine_type> > refineQueue;
	double integral_value;
	double move;
};

#endif /*INTEGRATION_HPP_*/
