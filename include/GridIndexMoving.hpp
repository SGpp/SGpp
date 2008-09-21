#ifndef GRIDINDEXMOVING_HPP_
#define GRIDINDEXMOVING_HPP_


#include "hash_map_config.hpp"

#include <iostream>
#include <sys/types.h>
#include <math.h>

template<class LT, class IT, int DIM>
class GridIndexMoving {
	
public:
	typedef LT level_type;
	typedef IT index_type;

	const static int dim = DIM;

    GridIndexMoving()
    {
    }

    GridIndexMoving(const GridIndexMoving<LT, IT, DIM>* o)
    {
        for(int d = 0; d < DIM; d++)
        {
            level[d] = o->level[d];
            index[d] = o->index[d];
            pos[d] = o->pos[d];
            base[d] = o->base[d];
        }
        rehash();
    }

	int getDim()
	{
		return DIM;
	}

	void set(int d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
		base[d] = pow(2.0, -static_cast<double>(l-1));
		pos[d] = i * base[d] * 0.5;
			
        rehash();
	}

	void push(int d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
		base[d] = pow(2.0, -static_cast<double>(l-1));
		pos[d] = i * base[d] * 0.5;

	}

	void get(int d, LT &l, IT &i) const
	{
		l = level[d];
		i = index[d];	
	}
	
	double getAbs(int d) const
	{
		return pos[d];	
	}

	double volume() const
	{
		double v = 1.0;
		for(int d = 0; d < DIM; d++)
		{
			v *= 0.5 * base[d];
		}
		return v;
	}
	
	double influence(int d, GridIndexMoving<LT, IT, DIM> &child)
	{
		return 0.5;
	}

	void move(int d, double p, double base)
	{
		this->pos[d] = p;
		this->base[d] = base;
	}

    void rehash()
    {
        size_t hash = 0xdeadbeef;
        for(int d = 0; d < DIM; d++)
        {
            //hash = level[d] + (hash << 6) + (hash << 16) - hash;
            //hash = index[d] + (hash << 6) + (hash << 16) - hash;
            //hash = (1<<level[d]) + index[d] + (hash << 6) + (hash << 16) - hash;
            hash = (1<<level[d]) + index[d] + hash*65599;
        }
        hash_value = hash;

/*
        size_t hash = 0xdeadbeef;
        char* p = (char*)(level);
        for(int d = 0; d < DIM*sizeof(LT); d++)
        {
        	hash = p[d] + (hash << 6) + (hash << 16) - hash;
        }

        p = (char*)(index);
        for(int d = 0; d < DIM*sizeof(IT); d++)
        {
        	hash = p[d] + (hash << 6) + (hash << 16) - hash;
        }
*/
        
    }

    size_t hash() const
    {
        return hash_value;
    }

    bool equals(GridIndexMoving<LT, IT, DIM> &rhs) const
    {
        for(int d = 0; d < DIM; d++)
        {
            if(level[d] != rhs.level[d])
            {
                return false;
            }
        }
        for(int d = 0; d < DIM; d++)
        {
            if(index[d] != rhs.index[d])
            {
                return false;
            }
        }
        return true;
    }

    GridIndexMoving<LT, IT, DIM>& operator= (const GridIndexMoving<LT, IT, DIM>& rhs)
    {
        for(int d = 0; d < DIM; d++)
        {
            level[d] = rhs.level[d];
            index[d] = rhs.index[d];
            base[d] = rhs.base[d];
            pos[d] = rhs.pos[d];
        }
        rehash();
        return *this;
    }

    void toString(std::ostream &stream)
    {
        stream << "[";
        for(int i = 0; i < DIM; i++)
        {
            if(i != 0)
            {
                stream << ",";
            }
            stream << " " << this->level[i];
            stream << ", " << this->index[i];
        }
        stream << " ]";
    }
    
    void printGrid(std::ostream &stream)
    {
        for(int i = 0; i < DIM; i++)
        {
			stream << pos[i] << " ";
        }
    }
    
    void print()
    {
    	
    	for(int i = 0; i < DIM; i++)
    	{
    		if(level[i] == 0)
    		{
    			std::cout << index[i];
    		}
    		else
    		{
    			std::cout << (pow(0.5, level[i])*index[i]);
    		}
    		std::cout << " ";
    	}
    }

private:
	LT level[DIM];
	IT index[DIM];
	double pos[DIM];
	double base[DIM];
    size_t hash_value;
	
};

/* no reference?
template<class LT, class IT, int DIM>
struct hash<GridIndex<LT, IT, DIM> > {
    size_t operator()(GridIndex<LT, IT, DIM> index) const {
        return index.hash();
    }
};
*/

template<class LT, class IT, int DIM>
struct hash<GridIndexMoving<LT, IT, DIM>* > {
    size_t operator()(GridIndexMoving<LT, IT, DIM>* index) const {
        return index->hash();
    }
};

/* no reference?
template<class LT, class IT, int DIM>
struct eqIndex<GridIndex<LT, IT, DIM> > {
    size_t operator()(GridIndex<LT, IT, DIM> &s1, GridIndex<LT, IT, DIM> &s2) const {
        return s1.equals(s2);
    }
};
*/

template<class LT, class IT, int DIM>
struct eqIndex<GridIndexMoving<LT, IT, DIM>* > {
    size_t operator()(GridIndexMoving<LT, IT, DIM>* s1, GridIndexMoving<LT, IT, DIM>* s2) const {
        return s1->equals(*s2);
    }
};


#endif /*GRIDINDEXMOVING_HPP_*/
