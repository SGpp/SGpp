// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <exception>
#include <list>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace sgpp {
namespace base {

HashGridStorage::HashGridStorage(size_t dimension)
    :  //  GridStorage(dim),
      dimension(dimension),
      list(),
      map(),
      algoDims(),
      boundingBox(new BoundingBox(dimension)),
      stretching(nullptr),
      bUseStretching(false) {
  for (size_t i = 0; i < dimension; i++) {
    algoDims.push_back(i);
  }
}

HashGridStorage::HashGridStorage(BoundingBox& creationBoundingBox)
    :  //  GridStorage(creationBoundingBox, creationBoundingBox.getDimensions()),
      dimension(creationBoundingBox.getDimension()),
      list(),
      map(),
      algoDims(),
      boundingBox(new BoundingBox(creationBoundingBox)),
      stretching(nullptr),
      bUseStretching(false) {
  // this look like a bug, creationBoundingBox not used
  for (size_t i = 0; i < dimension; i++) {
    algoDims.push_back(i);
  }
}

HashGridStorage::HashGridStorage(Stretching& creationStretching)
    :  //  : GridStorage(creationStretching, creationStretching.getDimensions()),
      dimension(creationStretching.getDimension()),
      list(),
      map(),
      algoDims(),
      boundingBox(nullptr),
      stretching(new Stretching(creationStretching)),
      bUseStretching(true) {
  // this look like a bug, creationBoundingBox not used
  for (size_t i = 0; i < dimension; i++) {
    algoDims.push_back(i);
  }
}

HashGridStorage::HashGridStorage(std::string& istr)
    :  //  : GridStorage(istr),
      dimension(0lu),
      list(),
      map(),
      algoDims() {
  std::istringstream istream;
  istream.str(istr);

  parseGridDescription(istream);

  for (size_t i = 0; i < dimension; i++) {
    algoDims.push_back(i);
  }
}

HashGridStorage::HashGridStorage(std::istream& istream)
    :  // GridStorage(istream),
      dimension(0lu),
      list(),
      map(),
      algoDims() {
  parseGridDescription(istream);

  for (size_t i = 0; i < dimension; i++) {
    algoDims.push_back(i);
  }
}

HashGridStorage::HashGridStorage(HashGridStorage& copyFrom)
    :  // GridStorage(copyFrom),
      dimension(copyFrom.dimension),
      list(),
      map(),
      algoDims(copyFrom.algoDims),
      boundingBox(copyFrom.bUseStretching ? nullptr : new BoundingBox(*copyFrom.boundingBox)),
      stretching(copyFrom.bUseStretching ? new Stretching(*copyFrom.stretching) : nullptr),
      bUseStretching(copyFrom.bUseStretching) {
  // copy gridpoints
  for (size_t i = 0; i < copyFrom.getSize(); i++) {
    this->insert(copyFrom[i]);
  }
}

void HashGridStorage::operator=(const HashGridStorage& other) {
  clear();

  if (bUseStretching) {
    delete stretching;
  } else {
    delete boundingBox;
  }

  dimension = other.dimension;
  algoDims = other.algoDims;
  bUseStretching = other.bUseStretching;

  if (other.bUseStretching) {
    stretching = new Stretching(*other.stretching);
  } else {
    boundingBox = new BoundingBox(*other.boundingBox);
  }

  for (size_t i = 0; i < other.getSize(); i++) {
    this->insert(other[i]);
  }
}

HashGridStorage::~HashGridStorage() {
  // delete all grid points
  if (bUseStretching) {
    delete stretching;
  } else {
    delete boundingBox;
  }

  for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
    delete *iter;
  }
}

void HashGridStorage::clear() {
  // delete all grid points
  for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
    delete *iter;
  }

  // remove all elements from hashmap
  map.clear();
  // remove all list entries
  list.clear();
}

std::vector<size_t> HashGridStorage::deletePoints(std::list<size_t>& removePoints) {
  point_pointer curPoint;
  std::vector<size_t> remainingPoints;
  size_t delCounter = 0;

  // sort list
  removePoints.sort();

  // DEBUG : print list points to delete, sorted
  // std::cout << std::endl << "List of points to delete, sorted" << std::endl;
  // for(std::list<size_t>::iterator iter = removePoints.begin();
  // iter != removePoints.end(); iter++)
  // {
  //   std::cout << " " << *iter << " ";
  // }
  // std::cout << std::endl;

  // Remove points with given indices for index vector and hashmap
  for (std::list<size_t>::iterator iter = removePoints.begin(); iter != removePoints.end();
       iter++) {
    size_t tmpIndex = *iter;
    size_t curPos = tmpIndex - delCounter;

    // GridPoint
    curPoint = list[curPos];

    // erase point
    delCounter++;
    map.erase(curPoint);
    list.erase(list.begin() + curPos);
  }

  // reset all entries in hash map and build list of remaining
  for (size_t i = 0; i < list.size(); i++) {
    curPoint = list[i];
    remainingPoints.push_back(map[curPoint]);
    map[curPoint] = i;
  }

  // reset the whole grid's leaf property in order
  // to guarantee a consistent grid
  recalcLeafProperty();

  // return indices of "surviver"
  return remainingPoints;
}

void HashGridStorage::unserializeNoAlgoDims(std::string& istr) {
  std::istringstream istream;
  istream.str(istr);

  parseGridDescription(istream);

  //    for (size_t i = 0; i < DIM; i++)
  //    {
  //      algoDims.push_back(i);
  //    }
}

std::string HashGridStorage::serialize(int version) const {
  std::ostringstream ostream;
  this->serialize(ostream, version);
  return ostream.str();
}

void HashGridStorage::serialize(std::ostream& ostream, int version) const {
  // Print version, dimensions and number of gridpoints
  ostream << version << " ";
  ostream << dimension << " ";
  ostream << list.size() << std::endl;

  // If BoundingBox used, write zero
  if (!bUseStretching) {
    ostream << std::scientific << 0 << std::endl;
    boundingBox->serialize(ostream, version);
  } else {
    // If analytic stretching, print the stretching type
    if (stretching->getStretchingMode() == "analytic") {
      ostream << std::scientific << 1 << std::endl;
    } else if (stretching->getStretchingMode() == "discrete") {
      // If discrete stretching, print the grid vector
      ostream << std::scientific << 2 << std::endl;
    }

    stretching->serialize(ostream, version);
  }

  // print the coordinates of the grid points
  for (grid_list_const_iterator iter = list.begin(); iter != list.end(); iter++) {
    (*iter)->serialize(ostream, version);
  }
}

std::string HashGridStorage::toString() const {
  std::ostringstream ostream;
  this->toString(ostream);
  return ostream.str();
}

void HashGridStorage::toString(std::ostream& stream) const {
  stream << "[";
  int i = 0;

  for (grid_map_const_iterator iter = map.begin(); iter != map.end(); iter++, i++) {
    if (i != 0) {
      stream << ",";
    }

    stream << " ";
    iter->first->toString(stream);
    stream << " -> " << iter->second;
  }

  stream << " ]";
}

size_t HashGridStorage::getSize() const { return map.size(); }

size_t HashGridStorage::getNumberOfInnerPoints() const {
  size_t innerPoints = 0;

  for (size_t p = 0; p < map.size(); p++) {
    if (list[p]->isInnerPoint()) innerPoints++;
  }

  return innerPoints;
}

size_t HashGridStorage::getDimension() const { return dimension; }

size_t HashGridStorage::insert(const point_type& index) {
  point_pointer insert = new HashGridPoint(index);
  list.push_back(insert);
  return (map[insert] = list.size() - 1);
}

void HashGridStorage::insert(point_type& index, std::vector<size_t>& insertedPoints) {
  index_t source_index;
  level_t source_level;

  if (!isContaining(index)) {
    // insert the current node
    size_t i = insert(index);
    insertedPoints.push_back(i);

    // insert all ancestors if they are missing
    for (size_t d = 0; d < dimension; d++) {
      // save level index in current dimension
      source_level = index.getLevel(d);
      source_index = index.getIndex(d);

      // go up to the parent node
      index.getParent(d);

      // insert all the parents until we find one which does already exist
      while (!isContaining(index)) {
        // insert the current node and all its missing ancestors
        insert(index, insertedPoints);
        index.getParent(d);
      }

      // reset index
      index.set(d, source_level, source_index);
    }
  }
}

void HashGridStorage::update(point_type& index, size_t pos) {
  if (pos < list.size()) {
    // Remove old element at pos
    point_pointer del = list[pos];
    map.erase(del);
    delete del;
    // Insert update
    point_pointer insert = new HashGridPoint(index);
    list[pos] = insert;
    map[insert] = pos;
  }
}

void HashGridStorage::deleteLast() {
  point_pointer del = list.back();
  map.erase(del);
  list.pop_back();
  delete del;
}

void HashGridStorage::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
  algoDims.clear();

  // throw an exception if there is
  if (newAlgoDims.size() > dimension) {
    throw generation_exception("There are more algorithmic dimensions than real dimensions!");
  }

  for (size_t i = 0; i < newAlgoDims.size(); i++) {
    algoDims.push_back(newAlgoDims[i]);
  }
}

void HashGridStorage::recalcLeafProperty() {
  point_pointer point;
  grid_map_iterator iter;
  size_t current_dim;
  point_type::level_type l;
  point_type::level_type i;
  bool isLeaf = true;

  // iterate through the grid
  for (iter = map.begin(); iter != map.end(); iter++) {
    point = iter->first;
    isLeaf = true;

    // iterate through the dimensions
    for (current_dim = 0; current_dim < dimension; current_dim++) {
      point->get(current_dim, l, i);

      if (l > 0) {
        // Test left child
        point->getLeftChild(current_dim);
        isLeaf = isLeaf && !isContaining(*point);

        // restore value for dimension
        point->set(current_dim, l, i);

        // Test right child
        point->getRightChild(current_dim);
        isLeaf = isLeaf && !isContaining(*point);
      } else {
        // Test level 0
        point->set(current_dim, 1, 1);
        isLeaf = isLeaf && !isContaining(*point);
      }

      // restore value for dimension
      point->set(current_dim, l, i);
    }

    point->setLeaf(isLeaf);
  }
}

// TODO(someone): this looks very fishy...
BoundingBox* HashGridStorage::getBoundingBox() {
  if (bUseStretching) {
    return stretching;
  } else {
    return boundingBox;
  }
}

Stretching* HashGridStorage::getStretching() { return stretching; }

void HashGridStorage::setBoundingBox(BoundingBox& boundingBox) {
  if (bUseStretching) {
    delete stretching;
  } else {
    delete this->boundingBox;
  }

  bUseStretching = false;
  this->boundingBox = new BoundingBox(boundingBox);
}

void HashGridStorage::setStretching(Stretching& stretching) {
  if (bUseStretching) {
    delete this->stretching;
  } else {
    delete boundingBox;
  }

  bUseStretching = true;
  this->stretching = new Stretching(stretching);
}

void HashGridStorage::getLevelIndexArraysForEval(DataMatrix& level, DataMatrix& index) {
  point_type::level_type curLevel;
  point_type::level_type curIndex;

  // Parallelization may lead to segfaults.... comment on your own risk
  //    #pragma omp parallel
  //    {
  //      #pragma omp for schedule (static) private(curLevel, curIndex)
  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);
      level.set(i, current_dim, static_cast<double>(1 << curLevel));
      index.set(i, current_dim, static_cast<double>(curIndex));
    }
  }

  //    }
}

void HashGridStorage::getLevelIndexArraysForEval(DataMatrixSP& level, DataMatrixSP& index) {
  point_type::level_type curLevel;
  point_type::level_type curIndex;

  // Parallelization may lead to segfaults.... comment on your own risk
  //    #pragma omp parallel
  //    {
  //      #pragma omp for schedule (static) private(curLevel, curIndex)
  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);
      level.set(i, current_dim, static_cast<float>(1 << curLevel));
      index.set(i, current_dim, static_cast<float>(curIndex));
    }
  }

  //    }
}

void HashGridStorage::getLevelForIntegral(DataMatrix& level) {
  point_type::level_type curLevel;
  point_type::level_type curIndex;

  // Parallelization may lead to segfaults.... comment on your own risk
  //    #pragma omp parallel
  //    {
  //      #pragma omp for schedule (static) private(curLevel, curIndex)
  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);
      level.set(i, current_dim, pow(2.0, static_cast<int>(-curLevel)));
    }
  }

  //    }
}

void HashGridStorage::getCoordinateArrays(DataMatrix& coordinates) {
  coordinates.resize(list.size(), dimension);

  base::DataVector x(dimension);
  for (size_t i = 0; i < list.size(); ++i) {
    getCoordinates(*list[i], x);
    coordinates.setRow(i, x);
  }
}

size_t HashGridStorage::getMaxLevel() const {
  point_type::level_type curLevel;
  point_type::level_type curIndex;
  point_type::level_type maxLevel;

  maxLevel = 0;

  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);

      if (curLevel > maxLevel) {
        maxLevel = curLevel;
      }
    }
  }

  return static_cast<size_t>(maxLevel);
}

void HashGridStorage::getLevelIndexMaskArraysForModEval(DataMatrix& level, DataMatrix& index,
                                                        DataMatrix& mask, DataMatrix& offset) {
  point_type::level_type curLevel;
  point_type::level_type curIndex;

  union IntMask {
    double d;
    uint64_t ui;
  } intMask;

  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);

      if (curLevel == 1) {
        level.set(i, current_dim, 0.0);
        index.set(i, current_dim, 0.0);
        intMask.ui = 0x0000000000000000;
        mask.set(i, current_dim, intMask.d);
        offset.set(i, current_dim, 1.0);
      } else if (curIndex == 1) {
        level.set(i, current_dim, (-1.0) * static_cast<double>(1 << curLevel));
        index.set(i, current_dim, 0.0);
        intMask.ui = 0x0000000000000000;
        mask.set(i, current_dim, intMask.d);
        offset.set(i, current_dim, 2.0);
      } else if (curIndex == static_cast<point_type::level_type>(((1 << curLevel) - 1))) {
        level.set(i, current_dim, static_cast<double>(1 << curLevel));
        index.set(i, current_dim, static_cast<double>(curIndex));
        intMask.ui = 0x0000000000000000;
        mask.set(i, current_dim, intMask.d);
        offset.set(i, current_dim, 1.0);
      } else {
        level.set(i, current_dim, static_cast<double>(1 << curLevel));
        index.set(i, current_dim, static_cast<double>(curIndex));
        intMask.ui = 0x8000000000000000;
        mask.set(i, current_dim, intMask.d);
        offset.set(i, current_dim, 1.0);
      }
    }
  }
}

void HashGridStorage::getLevelIndexMaskArraysForModEval(DataMatrixSP& level, DataMatrixSP& index,
                                                        DataMatrixSP& mask, DataMatrixSP& offset) {
  point_type::level_type curLevel;
  point_type::level_type curIndex;

  union IntMask {
    float f;
    uint32_t ui;
  } intMask;

  for (size_t i = 0; i < list.size(); i++) {
    for (size_t current_dim = 0; current_dim < dimension; current_dim++) {
      (list[i])->get(current_dim, curLevel, curIndex);

      if (curLevel == 1) {
        level.set(i, current_dim, 0.0);
        index.set(i, current_dim, 0.0);
        intMask.ui = 0x00000000;
        mask.set(i, current_dim, intMask.f);
        offset.set(i, current_dim, 1.0);
      } else if (curIndex == 1) {
        level.set(i, current_dim, (-1.0f) * static_cast<float>(1 << curLevel));
        index.set(i, current_dim, 0.0);
        intMask.ui = 0x00000000;
        mask.set(i, current_dim, intMask.f);
        offset.set(i, current_dim, 2.0);
      } else if (curIndex == static_cast<point_type::level_type>(((1 << curLevel) - 1))) {
        level.set(i, current_dim, static_cast<float>(1 << curLevel));
        index.set(i, current_dim, static_cast<float>(curIndex));
        intMask.ui = 0x00000000;
        mask.set(i, current_dim, intMask.f);
        offset.set(i, current_dim, 1.0);
      } else {
        level.set(i, current_dim, static_cast<float>(1 << curLevel));
        index.set(i, current_dim, static_cast<float>(curIndex));
        intMask.ui = 0x80000000;
        mask.set(i, current_dim, intMask.f);
        offset.set(i, current_dim, 1.0);
      }
    }
  }
}

void HashGridStorage::parseGridDescription(std::istream& istream) {
  int version;
  istream >> version;
  istream >> dimension;

  size_t num;
  istream >> num;

  // check whether grid was created with a version that is too new
  if (version > SERIALIZATION_VERSION) {
    if (version != 4) {
      std::ostringstream errstream;
      errstream << "Version of serialized grid (" << version
                << ") is too new. Max. recognized version is " << SERIALIZATION_VERSION << ".";
      throw generation_exception(errstream.str().c_str());
    }
  }

  // no bounding box, generate a trivial one
  if (version == 1 || version == 2) {
    // create a standard bounding box
    boundingBox = new BoundingBox(dimension);
  } else if (version == 3 || version == 4) {
    // read the bounding box
    // create a standard bounding box
    boundingBox = new BoundingBox(dimension);
    stretching = NULL;
    bUseStretching = false;
    BoundingBox1D tempBound;

    // reads the bounding box
    for (size_t i = 0; i < dimension; i++) {
      istream >> tempBound.leftBoundary;
      istream >> tempBound.rightBoundary;
      istream >> tempBound.bDirichletLeft;
      istream >> tempBound.bDirichletRight;

      boundingBox->setBoundary(i, tempBound);
    }
  } else if (version >= 5) {
    //      std::cout<<"Version 5 parse starts\n";
    int useStretching;
    BoundingBox1D tempBound;
    istream >> useStretching;

    if (useStretching == 0) {
      // BoundingBox
      boundingBox = new BoundingBox(dimension);
      stretching = nullptr;
      bUseStretching = false;
      boundingBox->unserialize(istream, version);
    } else {
      boundingBox = nullptr;
      stretching = new Stretching(dimension);
      bUseStretching = true;

      if (useStretching == 1) {
        // Stretching with analytic mode
        stretching->unserialize(istream, "analytic", version);
      } else if (useStretching == 2) {
        // Stretching with discrete Mode
        stretching->unserialize(istream, "discrete", version);
      } else {
        std::cout << "Unknown Container Id Given in parseGridDescription\n";
      }
    }
  }

  for (size_t i = 0; i < num; i++) {
    point_pointer index = new HashGridPoint(istream, version);
    list.push_back(index);
    map[index] = i;
  }

  // set's the grid point's leaf information which is not saved in version 1
  if (version == 1 || version == 4) {
    recalcLeafProperty();
  }
}

void HashGridStorage::getCoordinates(const HashGridPoint& point, DataVector& coordinates) const {
  coordinates.resize(dimension);

  for (size_t d = 0; d < dimension; d++) {
    coordinates[d] = getCoordinate(point, d);
  }
}

DataVector HashGridStorage::getCoordinates(const HashGridPoint& point) const {
  DataVector coordinates(dimension);

  for (size_t d = 0; d < dimension; d++) {
    coordinates[d] = getCoordinate(point, d);
  }

  return coordinates;
}

}  // namespace base
}  // namespace sgpp
