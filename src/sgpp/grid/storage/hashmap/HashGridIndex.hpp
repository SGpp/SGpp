/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHGRIDINDEX_HPP
#define HASHGRIDINDEX_HPP

#include "common/hash_map_config.hpp"
#include "data/DataVector.hpp"

#include "grid/common/BoundingBox.hpp"

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
        };

        if (version >= 2 && version != 4)
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
	 * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and rehashs the HashGridIndex object
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
	 * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and the Leaf property and rehashs the HashGridIndex object
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
	 * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and doesn't rehash the HashGridIndex object
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
	 * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and the Leaf property and doesn't rehash the HashGridIndex object
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
	 * gets level <i>l</i> and index <i>i</i> in dimension <i>d</i> by reference parameters
	 *
	 * @param d the dimension in which the ansatz function should be read
	 * @param l reference parameter for the level of the ansatz function
	 * @param i reference parameter for the index of the ansatz function
	 */
	void get(size_t d, LT &l, IT &i) const
	{
		l = level[d];
		i = index[d];
	}

	/**
	 * sets the Leaf option of this index
	 *
	 * @param isLeaf specifies if the current index is a leaf (i.e. has <b>no</b> child nodes) or not
	 */
	void setLeaf(bool isLeaf)
	{
		Leaf = isLeaf;
	}

	/**
	 * checks if this grid point has any successors in any dimension
	 *
	 * @return returns true if this grid point has no children otherwise false
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
	 * @todo (heinecke, should) rename to getCoord
	 */
	double abs(size_t d) const
	{
		return index[d] * pow(2.0, -static_cast<double>(level[d]));
	}

	/**
	 * determines the coordinate in a given dimension
	 *
	 * @param d the dimension in which the coordinate should be calculated
	 * @param q the intervals width in this dimension
	 * @param t the offset in this dimension
	 *
	 * @return the coordinate in the given dimension
	 */
	double getCoordBB(size_t d, double q, double t) const
	{
		return q*(index[d] * pow(2.0, -static_cast<double>(level[d]))) + t;
	}

	/**
	 * determines if the grid point is an inner grid point
	 *
	 * @return true if the grid point is an inner grid point
	 */
	bool isInnerPoint()
	{
		for (size_t d = 0; d < DIM; d++)
		{
			if (level[d] == 0)
				return false;
		}

		return true;
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
	 * when both f(x,y) and f(y,x) are false -> equalsSGLRBHash
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

#ifdef LARRABEENATIVE
    /**
     * checks whether this gridpoints is identical to another one
	 *
	 * Running under LARRABEE this is defined the way around, MSDN 2009:
	 * A binary predicate f(x,y) is a function object that has two
	 * argument objects x and y and a return value of true or false.
	 * An ordering imposed on a hash_map is a strict weak ordering
	 * if the binary predicate is irreflexive, antisymmetric,
	 * and transitive and if equivalence is transitive, where
	 * two objects x and y are defined to be equivalent
	 * when both f(x,y) and f(y,x) are false -> equalsSGLRBHash
     *
     * @param rhs reference the another HashGridIndex instance
     *
     * @return true if the gridpoints are identical otherwise false
     */
    bool equalsSGLRBHash(HashGridIndex<LT, IT> &rhs) const
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
	 * @todo who is "I"?????? --> (heinecke) I guess Joerg Blank ;-)
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
	 * operator to assign the current grid point with the values of another one
	 *
	 * @param rhs a reference to a HashGridIndex that contains the values that should be copied
	 *
	 * @todo (blank) generate working swig-wrapper
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
	 * Generates a string with level and index of the gridpoint.
     * The format is <tt>[l1, i1, l2, i2, ..., ld, id]</tt>.
     * Needed for Java compatibility.
	 *
	 * @returns string into which the gridpoint is written
	 */
	std::string toString()
	{
		std::ostringstream ostream;
		toString(ostream);

		return ostream.str();
	}

    /**
     * Generates a string with level and index of the gridpoint.
     * The format is <tt>[l1, i1, l2, i2, ..., ld, id]</tt>.
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
     * Sets the entries of DataVector p to the coordinates of the gridpoint with bounding box
     *
     * @param p the (result) DataVector p that should be overwritten
     * @param BB reference to BoundingBox Object, that stores all boundaries for all dimensions
     */
    void getCoordsBB(DataVector& p, BoundingBox& BB)
    {
		for(size_t i = 0; i < DIM; i++)
		{
			if(level[i] == 0)
			{
				p.set(i, (BB.getIntervalWidth(i)*index[i])+BB.getIntervalOffset(i));
			}
			else
			{
				p.set(i, (BB.getIntervalWidth(i)*(pow(0.5, static_cast<double>(level[i]))*index[i]))+BB.getIntervalOffset(i));
			}
		}
    }

    /**
     * Generates a string with all coordinates of the grid point.
     * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
     *
     * @return returns a string with the coordinates of the grid point separated by whitespace
     */
    std::string getCoordsString()
    {
    	std::stringstream return_stream;

    	// switch on scientific notation:
		//return_stream << std::scientific;

    	for(size_t i = 0; i < DIM; i++)
    	{
    		if(level[i] == 0)
    		{
    			return_stream << index[i];
    		}
    		else
    		{
    			return_stream << std::scientific << (pow(0.5, static_cast<double>(level[i]))*index[i]);
    		}

    		if (i < DIM-1)
    		{
    			return_stream << " ";
    		}
    	}

    	return return_stream.str();
    }

    /**
     * Generates a string with all coordinates of the grid point with bounding box
     * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
     *
     * This version scales the coordinates with q and t
     *
     * @param BB reference to BoundingBox Object, that stores all boundaries for all dimensionst
     *
     * @return returns a string with the coordinates of the grid point separated by whitespace
     */
    std::string getCoordsStringBB(BoundingBox& BB)
    {
    	std::stringstream return_stream;

    	for(size_t i = 0; i < DIM; i++)
    	{
    		if(level[i] == 0)
    		{
    			return_stream << std::scientific << (BB.getIntervalWidth(i)*index[i]) + BB.getIntervalOffset(i);
    		}
    		else
    		{
    			return_stream << std::scientific << (BB.getIntervalWidth(i)*(pow(0.5, static_cast<double>(level[i]))*index[i]))+BB.getIntervalOffset(i);
    		}

    		if (i < DIM-1)
    		{
    			return_stream << " ";
    		}
    	}

    	return return_stream.str();
    }

    /**
      * Returns the sum of the one-dimensional levels, i.e. @f$ |\vec{l}|_1 @f$.
      *
      * @return the sum of the one-dimensional levels
      */
     LT getLevelSum()
     {
         LT levelsum = 0;
         for (size_t i = 0; i < DIM; i++)
         {
             levelsum += level[i];
         }
         return levelsum;
     }

     /**
       * Returns the maximum of the one-dimensional levels, i.e. @f$ |\vec{l}|_\infty @f$.
       *
       * @return  the maximum of the one-dimensional levels
       */
      LT getLevelMax()
      {
          LT levelmax = 0;
          for (size_t i = 0; i < DIM; i++)
          {
              levelmax = max(levelmax, level[i]);
          }
          return levelmax;
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

#ifndef LARRABEENATIVE
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

#ifdef LARRABEENATIVE
#include <ext/hash_map>

/**
 * hash_compare function needed on Larrabee
 * for configuring the hash_hap
 */
template<class LT, class IT>
class LRBSGHasher<HashGridIndex<LT, IT>*> : public std::hash_compare<HashGridIndex<LT, IT>*>
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
		return s1->equalsSGLRBHash(*s2);
	}
};
#endif

}

#endif /* HASHGRIDINDEX_HPP */
