/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008-2009 Dirk Pflueger (dirk.pflueger@in.tum.de)           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                         */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef HASHGRIDINDEX_HPP
#define HASHGRIDINDEX_HPP

#include "common/hash_map_config.hpp"
#include "data/DataVector.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <cmath>

namespace sg
{

/**
 * This Class represents one Gridpoint.
 *
 * A Gridpoint is given by its
 * ansatzfunctions that are not zero in every dimension. Instances
 * of this class are members in the hashmap that represents the
 * whole grid.
 */
template<class LT, class IT>
class HashGridIndex {
public:
	typedef LT level_type;
	typedef IT index_type;

	/**
	 * Constructor of a n-Dim gridpoint
	 *
	 * @param dim the dimension of the gridpoint
	 */
    HashGridIndex(size_t dim) : DIM(dim), level(NULL), index(NULL)
    {
    	level = new level_type[dim];
    	index = new index_type[dim];
    	Leaf = false;
    }

    /**
     * Standard-Constructor
     */
    HashGridIndex() : DIM(0), level(NULL), index(NULL)
    {
    	Leaf = false;
    }

    /**
     * Copy-Constructor
     *
     * @param o constant pointer to HashGridIndex object
     */
    HashGridIndex(const HashGridIndex<LT, IT>* o) : DIM(o->DIM), level(NULL), index(NULL)
    {
    	level = new level_type[DIM];
    	index = new index_type[DIM];
    	Leaf = false;

        for(size_t d = 0; d < DIM; d++)
        {
            level[d] = o->level[d];
            index[d] = o->index[d];
        }
        Leaf = o->Leaf;
        rehash();
    }

    /**
     * Serialisation-Constructor
     *
     * @param istream instream object the contains the information about the gridpoint
     * @param version the serialization version of the file
     */
    HashGridIndex(std::istream& istream, int version) : DIM(0), level(NULL), index(NULL)
    {
    	size_t temp_leaf;

    	istream >> DIM;

    	level = new level_type[DIM];
    	index = new index_type[DIM];
    	Leaf = false;

        for(size_t d = 0; d < DIM; d++)
        {
			istream >> level[d];
			istream >> index[d];
        }

        if (version == 2)
        {
			// read leaf option
			istream >> temp_leaf;
			if (temp_leaf == 0)
			{
				Leaf = false;
			}
			else
			{
				Leaf = true;
			}
        }
        else
        {
        	Leaf = false;
        }

        rehash();
    }

    /**
     * Destructor
     */
	~HashGridIndex()
	{
		if(level)
		{
			delete [] level;
		}
		if(index)
		{
			delete [] index;
		}
	}

	/**
	 * Serialize this Gridpoint e.g. for a storage or checkpointing
	 *
	 * @param ostream outstream object to which the gridpoint's information is written
	 */
	void serialize(std::ostream& ostream)
	{
		ostream << DIM << std::endl;

        for(size_t d = 0; d < DIM; d++)
        {
			ostream << level[d] << " ";
			ostream << index[d] << " ";
        }
        ostream << std::endl;

        ostream << Leaf << std::endl;
	}

	/**
	 * Gets the dimension of the gridpoint
	 *
	 * @return the dimension of the gridpoint
	 */
	size_t dim() const
	{
		return DIM;
	}

	/**
	 * Sets the ansatzfunction in dimension <i>d</i> and rehashs the HashGridIndex object
	 *
	 * @param d the dimension in which the ansatzfunction is set
	 * @param l the level of the ansatzfunction
	 * @param i the index of the ansatzfunction
	 */
	void set(size_t d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
        rehash();
	}

	/**
	 * Sets the ansatzfunction in dimension <i>d</i> and the Leaf property and rehashs the HashGridIndex object
	 *
	 * @param d the dimension in which the ansatzfunction is set
	 * @param l the level of the ansatzfunction
	 * @param i the index of the ansatzfunction
	 * @param isLeaf specifies if this gridpoint has any childrens in any dimension
	 */
	void set(size_t d, LT l, IT i, bool isLeaf)
	{
		level[d] = l;
		index[d] = i;
		Leaf = isLeaf;
        rehash();
	}

	/**
	 * Sets the ansatzfunction in dimension <i>d</i> and doesn't rehash the HashGridIndex object
	 *
	 * @param d the dimension in which the ansatzfunction is set
	 * @param l the level of the ansatzfunction
	 * @param i the index of the ansatzfunction
	 */
	void push(size_t d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
	}

	/**
	 * Sets the ansatzfunction in dimension <i>d</i> and the Leaf property and doesn't rehash the HashGridIndex object
	 *
	 * @param d the dimension in which the ansatzfunction is set
	 * @param l the level of the ansatzfunction
	 * @param i the index of the ansatzfunction
	 * @param isLeaf specifies if this gridpoint has any childrens in any dimension
	 */
	void push(size_t d, LT l, IT i, bool isLeaf)
	{
		level[d] = l;
		index[d] = i;
		Leaf = isLeaf;
	}

	/**
	 * gets the ansatzfunction in given dimension by reference parameters
	 *
	 * @param d the dimension in which the ansatzfunction should be read
	 * @param l reference parameter for the level of the ansatzfunction
	 * @param i reference parameter for the index of the ansatzfunction
	 */
	void get(size_t d, LT &l, IT &i) const
	{
		l = level[d];
		i = index[d];
	}

	/**
	 * sets the Leaf option of this index
	 *
	 * @param isLeaf specifies if the current index is a leaf or not
	 */
	void setLeaf(bool isLeaf)
	{
		Leaf = isLeaf;
	}

	/**
	 * checks if this gridpoint has any succesors in any dimension
	 *
	 * @return returns true if this gridpoint has no children otherwise false
	 */
	bool isLeaf()
	{
		return Leaf;
	}

	/**
	 * determines the coordinate in a given dimension
	 *
	 * @param d the dimension in which the coordinate should be calculated
	 *
	 * @return the coordinate in the given dimension
	 * @todo rename to getCoord
	 */
	double abs(size_t d) const
	{
		return index[d] * pow(2.0, -static_cast<double>(level[d]));
	}

	/**
	 * gets a Pointer to the instance of the HashGridIndex Object
	 *
	 * @return Pointer to this instance
	 */
	HashGridIndex* getPointer()
	{
		return this;
	}

	/**
	 * rehashs the current gridpoint
	 */
    void rehash()
    {
        size_t hash = 0xdeadbeef;
        for(size_t d = 0; d < DIM; d++)
        {
            hash = (1<<level[d]) + index[d] + hash*65599;
        }
        hash_value = hash;
    }

    /**
     * gets the hash value of the current instance
     *
     * @return the hash value of the instance
     */
    size_t hash() const
    {
		return hash_value;
    }

    /**
     * checks whether this gridpoints is identical to another one
	 *
	 * Running under WINDOW this is defined the way around, MSDN 2009:
	 * A binary predicate f(x,y) is a function object that has two 
	 * argument objects x and y and a return value of true or false. 
	 * An ordering imposed on a hash_map is a strict weak ordering 
	 * if the binary predicate is irreflexive, antisymmetric, 
	 * and transitive and if equivalence is transitive, where 
	 * two objects x and y are defined to be equivalent 
	 * when both f(x,y) and f(y,x) are false -> equalsSGWinHash
     *
     * @param rhs reference the another HashGridIndex instance
     *
     * @return true if the gridpoints are identical otherwise false
     */
    bool equals(HashGridIndex<LT, IT> &rhs) const
    {
        for(size_t d = 0; d < DIM; d++)
        {
            if(level[d] != rhs.level[d])
            {
                return false;
            }
        }
        for(size_t d = 0; d < DIM; d++)
        {
            if(index[d] != rhs.index[d])
            {
                return false;
            }
        }
		return true;
    }

    /**
     * checks whether this gridpoints is identical to another one
	 *
	 * Running under WINDOW this is defined the way around, MSDN 2009:
	 * A binary predicate f(x,y) is a function object that has two 
	 * argument objects x and y and a return value of true or false. 
	 * An ordering imposed on a hash_map is a strict weak ordering 
	 * if the binary predicate is irreflexive, antisymmetric, 
	 * and transitive and if equivalence is transitive, where 
	 * two objects x and y are defined to be equivalent 
	 * when both f(x,y) and f(y,x) are false -> equalsSGWinHash
     *
     * @param rhs reference the another HashGridIndex instance
     *
     * @return true if the gridpoints are identical otherwise false
     */
#ifdef WINDOWS
    bool equalsSGWinHash(HashGridIndex<LT, IT> &rhs) const
    {
        for(size_t d = 0; d < DIM; d++)
        {
            if(level[d] != rhs.level[d])
            {
                return true;
            }
        }
        for(size_t d = 0; d < DIM; d++)
        {
            if(index[d] != rhs.index[d])
            {
                return true;
            }
        }
        return false;
    }
#endif

	/**
	 * This is just wrapper for operator= until I cant get swig to wrap it
	 *
	 * @param rhs a reference to a HashGridIndex that contains the values that should be copied
	 *
	 * @return returns a reference HashGridIndex
	 */
	HashGridIndex<LT, IT>& assign (const HashGridIndex<LT, IT>& rhs)
	{
		return this->operator=(rhs);
	}

	/**
	 * operator to assign the current gridpoint with the values of another one
	 *
	 * @param rhs a reference to a HashGridIndex that contains the values that should be copied
	 *
	 * @todo: generate working swig-wrapper
	 *
	 * @return returns a reference HashGridIndex
	 */
    HashGridIndex<LT, IT>& operator= (const HashGridIndex<LT, IT>& rhs)
    {
    	if(this == &rhs)
    	{
    		return *this;
    	}

    	if(DIM != rhs.DIM)
    	{

    		if(level)
    		{
    			delete [] level;
    		}
    		if(index)
    		{
    			delete [] index;
    		}

    		DIM = rhs.DIM;

    		level = new level_type[DIM];
    		index = new index_type[DIM];
    	}

        for(size_t d = 0; d < DIM; d++)
        {
            level[d] = rhs.level[d];
            index[d] = rhs.index[d];
        }
        Leaf = rhs.Leaf;

        rehash();
        return *this;
    }

    /**
     * generates a string with ansatzfunctions of the gridpoint
     *
     * @param stream reference to a output stream
     */
    void toString(std::ostream &stream)
    {
        stream << "[";
        for(size_t i = 0; i < DIM; i++)
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

    /**
     * Sets the entries of DataVector p to the coordinates of the gridpoint
     *
     * @param p the (result) DataVector p that should be overwritten
     */
    void getCoords(DataVector &p)
    {
      for(size_t i = 0; i < DIM; i++)
    	{
	  if(level[i] == 0)
	    {
	      p.set(i, index[i]);
	    }
	  else
	    {
	      p.set(i, pow(0.5, static_cast<double>(level[i]))*index[i]);
	    }
	}      
    }


    /**
     * generates a string with all coordinates of the gridpoint
     *
     * @return returns a string with the coordinates of the gridpoint separated by whitespace
     * @todo rename to getCoordString()
     */
    std::string getCoordinates()
    {
    	std::stringstream return_stream;

    	for(size_t i = 0; i < DIM; i++)
    	{
    		if(level[i] == 0)
    		{
    			return_stream << index[i];
    		}
    		else
    		{
    			return_stream << (pow(0.5, static_cast<double>(level[i]))*index[i]);
    		}

    		if (i < DIM-1)
    		{
    			return_stream << " ";
    		}
    	}

    	return return_stream.str();
    }

private:
	/// the dimension of the gridpoint
	size_t DIM;
	/// pointer to array that stores the ansatzfunctions' level
	LT* level;
	/// pointer to array that stores the ansatzfunctions' indecies
	IT* index;
	/// stores if this gridpoint is a leaf
	bool Leaf;
	/// stores the hashvalue of the gridpoint
	size_t hash_value;
};

#ifndef WINDOWS
template<class LT, class IT>
struct hash<HashGridIndex<LT, IT>* > {
    size_t operator()(HashGridIndex<LT, IT>* index) const {
        return index->hash();
    }
};

template<class LT, class IT>
struct eqIndex<HashGridIndex<LT, IT>* > {
    size_t operator()(HashGridIndex<LT, IT>* s1, HashGridIndex<LT, IT>* s2) const {
        return s1->equals(*s2);
    }
};
#endif
#ifdef WINDOWS
#include <hash_map>

/**
 * hash_compare function needed on Windows boxes
 * for configuring the hash_hap
 */
template<class LT, class IT>
class WinSGHasher<HashGridIndex<LT, IT>*> : public stdext::hash_compare<HashGridIndex<LT, IT>*>
{
public:
	/**
	 * operator that calculates the hash values
	 *
	 * @param index pointer to an element of the hash_map
	 * @return the hash value
	 */
	size_t operator() (HashGridIndex<LT, IT>* index) const
	{
		return index->hash();
	}

	/**
	 * operator that compares to elements in the hash_map
	 *
	 * @param s1 first element
	 * @param s2 second element
	 * @return returns true if s1 == s2, otherwise false
	 */
	bool operator() (HashGridIndex<LT, IT>* s1, HashGridIndex<LT, IT>* s2) const
	{
		return s1->equalsSGWinHash(*s2);
	}
};
#endif
}

#endif /* HASHGRIDINDEX_HPP */
