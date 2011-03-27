/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHGRIDSTORAGE_HPP
#define HASHGRIDSTORAGE_HPP

#include "common/hash_map_config.hpp"

#include "exception/generation_exception.hpp"

#include "grid/storage/hashmap/HashGridIndex.hpp"
#include "grid/storage/hashmap/HashGridIterator.hpp"
#include "grid/storage/hashmap/SerializationVersion.hpp"

#include "grid/common/BoundingBox.hpp"

#include <memory>
#include <string>
#include <sstream>
#include <exception>

namespace sg
{

/**
 * Generic hash table based index storage.
 */
template<typename GIT>
class HashGridStorage
{
public:
    typedef GIT index_type;
    typedef GIT* index_pointer;
#ifndef USETRONE
#ifndef LARRABEENATIVE
    typedef std::hash_map<index_pointer, size_t, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
#endif
#ifdef LARRABEENATIVE
	typedef std::hash_map<index_pointer, size_t, LRBSGHasher<index_pointer> > grid_map;
#endif
#endif
#ifdef USETRONE
    typedef std::tr1::unordered_map<index_pointer, size_t, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
#endif
    typedef typename grid_map::iterator grid_map_iterator;
    typedef typename grid_map::const_iterator grid_map_const_iterator;

	typedef std::vector<index_pointer> grid_list;
	typedef typename grid_list::iterator grid_list_iterator;

	typedef HashGridIterator<GIT> grid_iterator;

	/**
	 * Constructor
	 *
	 * initializes the boundingBox with a trivial cube
	 *
	 * @param dim the dimension of the sparse grid
	 */
	HashGridStorage(size_t dim) : DIM(dim), list(), map()
	{
		boundingBox = new BoundingBox(DIM);
	}

	/**
	 * Constructor
	 *
	 * initializes the boundingBox with a reference to another boundingbox
	 *
	 * @param creationBoundingBox reference to bounding box object that describes the grid boundaries
	 */
	HashGridStorage(BoundingBox& creationBoundingBox) : DIM(), list(), map()
	{
		boundingBox = new BoundingBox(creationBoundingBox);
		DIM = boundingBox->getDimensions();
	}

	/**
	 * Constructor that reads the data from a string
	 *
	 * @param istr the string that contains the data
	 */
	HashGridStorage(std::string& istr) : DIM(0), list(), map()
	{
    	std::istringstream istream;
    	istream.str(istr);

		parseGridDescription(istream);
	}

	/**
	 * Constructor that reads the data from an input stream
	 *
	 * @param istream the inputstream that contains the data
	 */
	HashGridStorage(std::istream& istream) : DIM(0), list(), map()
	{
		parseGridDescription(istream);
	}


	/**
	 * Desctructor
	 */
	~HashGridStorage()
	{
		// delete all grid points
		for(grid_list_iterator iter = list.begin(); iter != list.end(); iter++)
		{
			delete *iter;
		}

		// delete hash map of grid indices
		// --> is auto deleted, because it's a class's element

		// delete the grid's bounding box
		delete boundingBox;
	}

	/**
	 * deletes all grid points in the storage
	 */
	void emptyStorage()
	{
		// delete all grid points
		for(grid_list_iterator iter = list.begin(); iter != list.end(); iter++)
		{
			delete *iter;
		}

		// remove all elements from hashmap
		map.clear();
		// remove all list entries
		list.clear();
	}

	/**
	 * serialize the gridstorage into a string
	 *
	 * @return a string the contains all gridstorage information
	 */
	std::string serialize()
	{
		std::ostringstream ostream;
		this->serialize(ostream);
		return ostream.str();
	}

	/**
	 * serialize the gridstorage into a stream
	 *
	 * @param ostream reference to a stream into that all gridstorage information is written
	 */
	void serialize(std::ostream& ostream)
	{
		DimensionBoundary tempBound;

		// Print version, dimensions and number of gridpoints
		ostream << SERIALIZATION_VERSION << " ";
		ostream << DIM << " ";
		ostream << list.size() << std::endl;

		// Print the bounding box
		for (size_t i = 0; i < DIM; i++)
		{
			tempBound = boundingBox->getBoundary(i);
			ostream << std::scientific << tempBound.leftBoundary << " " << tempBound.rightBoundary << " " << tempBound.bDirichletLeft << " " << tempBound.bDirichletRight << " ";
		}
		ostream << std::endl;

		// print the coordinates of the grid points
		for(grid_list_iterator iter = list.begin(); iter != list.end(); iter++)
		{
			(*iter)->serialize(ostream);
		}
	}

    /**
	 * serialize the gridstorage's gridpoints into a stream
	 *
	 * @return returns the string that contains all gridpoint information
	 */
    std::string toString()
	{
		std::ostringstream ostream;
		this->toString(ostream);
		return ostream.str();
	}

	/**
	 * serialize the gridstorage's gridpoints into a stream
	 *
	 * @param stream reference to a stream into that all gridpoint information is written
	 */
    void toString(std::ostream& stream)
    {
        stream << "[";
        int i = 0;
       	grid_map_iterator iter;
        for(iter = map.begin(); iter != map.end(); iter++, i++)
        {
            if(i != 0)
            {
                stream << ",";
            }
            stream << " ";
            iter->first->toString(stream);
            stream << " -> " << iter->second;
        }

        stream << " ]";
    }

    /**
     * gets the size of the hashmap
     *
     * @return returns the size of the hashmap
     */
    size_t size() const
    {
        return map.size();
    }

    /**
     * gets the number of inner grid points
     *
     * @return the number of inner grid points
     */
    size_t getNumInnerPoints()
    {
    	size_t innerPoints = 0;

    	for (size_t p = 0; p < map.size(); p++)
    	{
    		if (list[p]->isInnerPoint())
    			innerPoints++;
    	}

    	return innerPoints;
    }

    /**
     * gets the dimension of the grid
     *
     * @return the dimension of the grid stored in this GridStorage object
     */
	size_t dim() const
	{
		return DIM;
	}

	/**
	 * gets the sequence number for given gridpoint by its index
	 *
	 * @param index gridindex object
	 *
	 * @return sequence number of the given index
	 */
	size_t& operator[](index_pointer index)
	{
		return map[index];
	}

	/**
	 * gets the index number for given gridpoint by its sequence number
	 *
	 * @param seq the sequence number of the index
	 *
	 * @return gridindex object (reference)
	 */
	index_pointer& operator[](size_t seq)
	{
		return list[seq];
	}

	/**
	 * gets the index number for given gridpoint by its sequence number
	 *
	 * @param seq the sequence number of the index
	 *
	 * @return gridindex object (pointer)
	 */
	GIT* get(size_t seq)
	{
		return list[seq];
	}

	/**
	 * insert a new index into map
	 *
	 * @param index reference to the index that should be inserted
	 *
	 * @return
	 */
	size_t insert(index_type &index)
	{
		index_pointer insert = new GIT(&index);
		list.push_back(insert);
		return (map[insert] = this->seq() - 1);
	}

	/**
	 * creates a pointer to index from a reference to index by creating
	 * a new instance of an index object
	 *
	 * @param index address of index object
	 *
	 * @return pointer to new index object
	 */
	index_pointer create(index_type &index)
	{
		index_pointer insert = new GIT(&index);
		return insert;
	}

	/**
	 * removes an index from gridstorage
	 *
	 * @param index pointer to index that should be removed
	 */
	void destroy(index_pointer index)
	{
		delete index;
	}

	/**
	 * stores a given index in the hashmap
	 *
	 * @param index pointer to index that should be stored
	 *
	 * @return sequence number
	 */
	unsigned int store(index_pointer index)
	{
		list.push_back(index);
		return (map[index] = this->seq() - 1);
	}

	/**
	 * sets the iterator to a given index
	 *
	 * @param index the index to which the cursor should be moved
	 */
	grid_map_iterator find(index_pointer index)
	{
		return map.find(index);
	}

	/**
	 * set iterator to the first position in the map
	 */
	grid_map_iterator begin()
	{
		return map.begin();
	}

	/**
	 * sets the iterator to last position in the map
	 */
	grid_map_iterator end()
	{
		return map.end();
	}

	/**
	 * Tests if index is in the storage
	 *
	 * @param index pointer to index that should be tested
	 *
	 * @return true if the index is in the storage
	 */
	bool has_key(GIT* index)
	{
		return map.find(index) != map.end();
	}

	/**
	 * Sets the index to the left level zero parent
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void left_levelzero(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, 0, 0);
	}

	/**
	 * Sets the index to the right level zero parent
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void right_levelzero(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, 0, 1);
	}

	/**
	 * Sets the index to the left child
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void left_child(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, l + 1, 2 * i - 1);
	}

	/**
	 * Sets the index to the right child
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void right_child(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, l + 1, 2 * i + 1);
	}

	/**
	 * Resets the index to the top level in direction d
	 *
	 * @param index pointer to index the should be modified
	 * @param d the dimension in which the modification is taken place
	 */
	void top(GIT* index, size_t d)
	{
		index->set(d, 1, 1);
	}

	/**
	 * Gets the seq number for index
	 *
	 * @param index pointer to index which sequence number should be determined
	 *
	 * @return the seq number for index
	 */
	size_t seq(GIT* index)
	{
		grid_map_iterator iter = map.find(index);
		if(iter != map.end())
		{
			return iter->second;
		}
		else
		{
			return map.size() + 1;
		}
	}

	/**
	 * Tests if seq number does not point to a valid grid index
	 *
	 * @param s sequence number that should be tested
	 *
	 * @return true if we are not EOF
	 */
	bool end(size_t s)
	{
		return s > map.size();
	}

	/**
	 * Recalculates the leaf-property of every grid point.
	 * This might be useful in case of a grid unserialization
	 */
	void recalcLeafProperty()
	{
       	index_pointer point;
		grid_map_iterator iter;
		size_t current_dim;
		typename index_type::level_type l;
		typename index_type::level_type i;
		bool isLeaf = true;

		// iterate through the grid
		for(iter = map.begin(); iter != map.end(); iter++)
        {
			point = iter->first;
			isLeaf = true;

			// iterate through the dimensions
			for (current_dim = 0; current_dim < DIM; current_dim++)
			{
				point->get(current_dim, l, i);

				if (l > 0)
				{
					// Test left child
					left_child(point, current_dim);
					isLeaf = isLeaf && !has_key(point);

					// restore value for dimension
					point->set(current_dim, l, i);

					// Test right child
					right_child(point, current_dim);
					isLeaf = isLeaf && !has_key(point);
				}
				else
				{
					// Test level 0
					point->set(current_dim, 1, 1);
					isLeaf = isLeaf && !has_key(point);
				}

				// restore value for dimension
				point->set(current_dim, l, i);
			}

			point->setLeaf(isLeaf);
        }
	}

	/**
	 * get the bounding box of the current grid
	 *
	 * @return returns a pointer to GridStorage's bounding box
	 */
	BoundingBox* getBoundingBox()
	{
		return boundingBox;
	}

	/**
	 * sets the bounding box of the current grid
	 *
	 * @param bb bounding box to which the GridStorage's pointer is set
	 */
	void setBoundingBox(BoundingBox& bb)
	{
		delete boundingBox;
		boundingBox = new BoundingBox(bb);
	}

protected:
	/**
	 * returns the next sequence numbers
	 *
	 * @return returns the next sequence numbers
	 */
    size_t seq()
    {
        return list.size();
    }

private:
	/// the dimension of the grid
	size_t DIM;
	/// the gridpoints
	grid_list list;
	/// the indices of the grid points
    grid_map map;
    /// the grids bounding box
    BoundingBox* boundingBox;

    /**
     * Parses the gird's information (grid points, dimensions, bounding box) from a string stream
     *
     * @param istream the string stream that contains the information
     */
    void parseGridDescription(std::istream& istream)
    {
    	int version;
    	istream >> version;

    	istream >> DIM;

    	// create a standard bounding box
    	boundingBox = new BoundingBox(DIM);

    	size_t num;
    	istream >> num;

		// check whether grid was created with a version that is too new
		if (version > SERIALIZATION_VERSION)
		{
			if (version != 4)
			{
				//	    throw generation_exception("Version of serialized grid is too new. Max. recognized version is "+string(SERIALIZATION_VERSION));
				std::ostringstream errstream;
				errstream << "Version of serialized grid (" << version << ") is too new. Max. recognized version is " << SERIALIZATION_VERSION << ".";
				throw generation_exception(errstream.str().c_str());
			}
		}

    	// read the bounding box
    	if (version == 3 || version == 4)
    	{
    		DimensionBoundary tempBound;

    		// reads the bounding box
    		for (size_t i = 0; i < DIM; i++)
    		{
    			istream >> tempBound.leftBoundary;
    			istream >> tempBound.rightBoundary;
    			istream >> tempBound.bDirichletLeft;
    			istream >> tempBound.bDirichletRight;

    			boundingBox->setBoundary(i, tempBound);
    		}
    	}

    	for(size_t i = 0; i < num; i++)
    	{
    		index_pointer index = new GIT(istream, version);
    		list.push_back(index);
    		map[index] = i;
    	}

    	// set's the grid point's leaf information which is not saved in version 1
    	if (version == 1 || version == 4)
    	{
    		recalcLeafProperty();
    	}
    }
};

}

#endif /* HASHGRIDSTORAGE_HPP */
