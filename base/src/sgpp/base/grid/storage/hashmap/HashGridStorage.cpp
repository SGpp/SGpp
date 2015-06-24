
#include <memory>
#include <exception>
#include <typeinfo>
#include <unordered_map>

#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>


#include <sgpp/base/exception/generation_exception.hpp>


namespace SGPP {
  namespace base {

    HashGridStorage::HashGridStorage(size_t dim) :
      //  GridStorage(dim),
      DIM(dim)
      , list()
      , map()
      , algoDims()
      , boundingBox(new BoundingBox(dim))
      , stretching(nullptr)
      , bUseStretching(false) {
      for (size_t i = 0; i < DIM; i++) {
        algoDims.push_back(i);
      }
    }

    HashGridStorage::HashGridStorage(BoundingBox& creationBoundingBox) :
      //  GridStorage(creationBoundingBox, creationBoundingBox.getDimensions()),
      DIM(creationBoundingBox.getDimensions())
      , list()
      , map()
      , algoDims()
      , boundingBox(new BoundingBox(creationBoundingBox))
      , stretching(nullptr)
      , bUseStretching(false) {
      //TODO gerrit: this look like a bug, creationBoundingBox not used
      for (size_t i = 0; i < DIM; i++) {
        algoDims.push_back(i);
      }
    }

    HashGridStorage::HashGridStorage(Stretching& creationStretching) :
      //  : GridStorage(creationStretching, creationStretching.getDimensions()),
      DIM(creationStretching.getDimensions())
      , list()
      , map()
      , algoDims()
      , boundingBox(nullptr)
      , stretching(new Stretching(creationStretching))
      , bUseStretching(true) {
      //TODO gerrit: this look like a bug, creationBoundingBox not used
      for (size_t i = 0; i < DIM; i++) {
        algoDims.push_back(i);
      }
    }

    HashGridStorage::HashGridStorage(std::string& istr) :
      //  : GridStorage(istr),
      DIM(0lu)
      , list()
      , map()
      , algoDims() {
      std::istringstream istream;
      istream.str(istr);

      parseGridDescription(istream);

      for (size_t i = 0; i < DIM; i++) {
        algoDims.push_back(i);
      }
    }


    HashGridStorage::HashGridStorage(std::istream& istream) :
      // GridStorage(istream),
      DIM(0lu)
      , list()
      , map()
      , algoDims() {
      parseGridDescription(istream);

      for (size_t i = 0; i < DIM; i++) {
        algoDims.push_back(i);
      }
    }

    HashGridStorage::HashGridStorage(HashGridStorage& copyFrom) :
      // GridStorage(copyFrom),
      DIM(copyFrom.DIM)
      , list()
      , map()
      , algoDims(copyFrom.algoDims)
      , boundingBox(copyFrom.bUseStretching ? nullptr : new BoundingBox(*copyFrom.boundingBox))
      , stretching(copyFrom.bUseStretching ? new Stretching(*copyFrom.stretching) : nullptr)
      , bUseStretching(copyFrom.bUseStretching) {
      // copy gridpoints
      for (size_t i = 0; i < copyFrom.size(); i++) {
        this->insert(*(copyFrom[i]));
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

    void
    HashGridStorage::emptyStorage() {
      // delete all grid points
      for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
        delete *iter;
      }

      // remove all elements from hashmap
      map.clear();
      // remove all list entries
      list.clear();
    }

    std::vector<size_t>
    HashGridStorage::deletePoints(std::list<size_t>& removePoints) {
      index_pointer curPoint;
      std::vector<size_t> remainingPoints;
      size_t delCounter = 0;

      // sort list
      removePoints.sort();

      //DEBUG : print list points to delete, sorted
      //std::cout << std::endl << "List of points to delete, sorted" << std::endl;
      //for(std::list<size_t>::iterator iter = removePoints.begin(); iter != removePoints.end(); iter++)
      //{
      //  std::cout << " " << *iter << " ";
      //}
      //std::cout << std::endl;

      // Remove points with given indices for index vector and hashmap
      for (std::list<size_t>::iterator iter = removePoints.begin();
           iter != removePoints.end(); iter++) {
        size_t tmpIndex = *iter;
        size_t curPos = tmpIndex - delCounter;

        // GridIndex
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

    void
    HashGridStorage::unserialize_noAlgoDims(std::string& istr) {
      std::istringstream istream;
      istream.str(istr);

      parseGridDescription(istream);

      //    for (size_t i = 0; i < DIM; i++)
      //    {
      //      algoDims.push_back(i);
      //    }
    }

    std::string
    HashGridStorage::serialize() {
      std::ostringstream ostream;
      this->serialize(ostream);
      return ostream.str();
    }

    void
    HashGridStorage::serialize(std::ostream& ostream) {
      DimensionBoundary tempBound;

      // Print version, dimensions and number of gridpoints
      ostream << SERIALIZATION_VERSION << " ";
      ostream << DIM << " ";
      ostream << list.size() << std::endl;

      //If BoundingBox used, write zero
      if (!bUseStretching) {
        ostream << std::scientific << 0 << std::endl;

        // Print the bounding box
        for (size_t i = 0; i < DIM; i++) {
          tempBound = boundingBox->getBoundary(i);
          ostream << std::scientific << tempBound.leftBoundary << " "
                  << tempBound.rightBoundary << " " << tempBound.bDirichletLeft << " "
                  << tempBound.bDirichletRight << " ";
        }

        ostream << std::endl;
      } else {
        //If analytic stretching, print the stretching type
        if (*(stretching->getStretchingMode()) == "analytic") {
          ostream << std::scientific << 1 << std::endl;

          // Print the bounding box
          for (size_t i = 0; i < DIM; i++) {
            tempBound = stretching->getBoundary(i);
            ostream << std::scientific << tempBound.leftBoundary << " "
                    << tempBound.rightBoundary << " " << tempBound.bDirichletLeft
                    << " " << tempBound.bDirichletRight << " ";
          }

          ostream << std::endl;
          Stretching1D str1d;
          int stretchingType = 0;

          /*
           * Write stretching type if
           * id: 1
           * log: 2
           * sinh:3
           */
          for (size_t i = 0; i < DIM; i++) {
            str1d = stretching->getStretching1D(i);

            if (str1d.type == "id") {
              stretchingType = 1;
            } else if (str1d.type == "log") {
              stretchingType = 2;
            } else if (str1d.type == "sinh") {
              stretchingType = 3;
            }

            ostream << std::scientific << stretchingType << " " << str1d.x_0
                    << " " << str1d.xsi << std::endl;

          }
        }
        //If discrete stretching, print the grid vector
        else if (*(stretching->getStretchingMode()) == "discrete") {
          ostream << std::scientific << 2 << std::endl;

          // Print the bounding box
          for (size_t i = 0; i < DIM; i++) {
            tempBound = stretching->getBoundary(i);
            ostream << std::scientific << tempBound.leftBoundary << " "
                    << tempBound.rightBoundary << " " << tempBound.bDirichletLeft
                    << " " << tempBound.bDirichletRight << " ";
          }

          ostream << std::endl;
          std::vector<float_t>* vec = stretching->getDiscreteVector(true);
          int* vecLevel = stretching->getDiscreteVectorLevel();

          for (size_t i = 0; i < DIM; i++) {
            ostream << std::scientific << vecLevel[i] << std::endl;

            for (size_t j = 0; j < vec[i].size(); j++) {
              ostream << std::scientific << vec[i][j] << " ";
            }

            ostream << std::endl;
          }
        }
      }

      // print the coordinates of the grid points
      for (grid_list_iterator iter = list.begin(); iter != list.end(); iter++) {
        (*iter)->serialize(ostream);
      }
    }

    std::string
    HashGridStorage::toString() {
      std::ostringstream ostream;
      this->toString(ostream);
      return ostream.str();
    }

    void
    HashGridStorage::toString(std::ostream& stream) {
      stream << "[";
      int i = 0;
      grid_map_iterator iter;

      for (iter = map.begin(); iter != map.end(); iter++, i++) {
        if (i != 0) {
          stream << ",";
        }

        stream << " ";
        iter->first->toString(stream);
        stream << " -> " << iter->second;
      }

      stream << " ]";
    }

    size_t
    HashGridStorage::size() const {
      return map.size();
    }

    size_t
    HashGridStorage::getNumInnerPoints() const {
      size_t innerPoints = 0;

      for (size_t p = 0; p < map.size(); p++) {
        if (list[p]->isInnerPoint())
          innerPoints++;
      }

      return innerPoints;
    }

    size_t
    HashGridStorage::dim() const {
      return DIM;
    }

    size_t
    HashGridStorage::insert(index_type& index) {
      index_pointer insert = new HashGridIndex(&index);
      list.push_back(insert);
      return (map[insert] = this->seq() - 1);
    }

    void
    HashGridStorage::update(index_type& index, size_t pos) {
      if (pos < seq()) {
        // Remove old element at pos
        index_pointer del = list[pos];
        map.erase(del);
        delete del;
        // Insert update
        index_pointer insert = new HashGridIndex(&index);
        list[pos] = insert;
        map[insert] = pos;
      }
    }

    void
    HashGridStorage::deleteLast() {
      index_pointer del = list.back();
      map.erase(del);
      list.pop_back();
      delete del;
    }


    void
    HashGridStorage::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
      algoDims.clear();

      // throw an exception if there is
      if (newAlgoDims.size() > DIM) {
        throw generation_exception(
          "There are more algorithmic dimensions than real dimensions!");
      }

      for (size_t i = 0; i < newAlgoDims.size(); i++) {
        algoDims.push_back(newAlgoDims[i]);
      }
    }

    void
    HashGridStorage::recalcLeafProperty() {
      index_pointer point;
      grid_map_iterator iter;
      size_t current_dim;
      index_type::level_type l;
      index_type::level_type i;
      bool isLeaf = true;

      // iterate through the grid
      for (iter = map.begin(); iter != map.end(); iter++) {
        point = iter->first;
        isLeaf = true;

        // iterate through the dimensions
        for (current_dim = 0; current_dim < DIM; current_dim++) {
          point->get(current_dim, l, i);

          if (l > 0) {
            // Test left child
            left_child(point, current_dim);
            isLeaf = isLeaf && !has_key(point);

            // restore value for dimension
            point->set(current_dim, l, i);

            // Test right child
            right_child(point, current_dim);
            isLeaf = isLeaf && !has_key(point);
          } else {
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

    BoundingBox*
    HashGridStorage::getBoundingBox() {
      return boundingBox;
    }

    Stretching*
    HashGridStorage::getStretching() {
      return stretching;
    }

    void
    HashGridStorage::setBoundingBox(BoundingBox& bb) {
      if (!bUseStretching) {
        delete boundingBox;
      }

      if (bUseStretching) {
        delete stretching;
      }

      boundingBox = new BoundingBox(bb);
      bUseStretching = false;
    }

    void
    HashGridStorage::setStretching(Stretching& bb) {
      if (bUseStretching) {
        delete stretching;
      }

      if (!bUseStretching) {
        delete boundingBox;
      }

      bUseStretching = true;
      stretching = new Stretching(bb);
    }

    void
    HashGridStorage::getLevelIndexArraysForEval(DataMatrix& level, DataMatrix& index) {
      index_type::level_type curLevel;
      index_type::level_type curIndex;

      // Parallelization may lead to segfaults.... comment on your own risk
      //    #pragma omp parallel
      //    {
      //      #pragma omp for schedule (static) private(curLevel, curIndex)
      for (size_t i = 0; i < list.size(); i++) {
        for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
          (list[i])->get(current_dim, curLevel, curIndex);
          level.set(i, current_dim, static_cast<float_t>(1 << curLevel));
          index.set(i, current_dim, static_cast<float_t>(curIndex));
        }
      }

      //    }
    }

    void
    HashGridStorage::getLevelForIntegral(DataMatrix& level) {
      index_type::level_type curLevel;
      index_type::level_type curIndex;

      // Parallelization may lead to segfaults.... comment on your own risk
      //    #pragma omp parallel
      //    {
      //      #pragma omp for schedule (static) private(curLevel, curIndex)
      for (size_t i = 0; i < list.size(); i++) {
        for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
          (list[i])->get(current_dim, curLevel, curIndex);
          level.set(i, current_dim, pow(2.0, static_cast<int>(-curLevel)));
        }
      }

      //    }
    }


    size_t
    HashGridStorage::getMaxLevel() const {
      index_type::level_type curLevel;
      index_type::level_type curIndex;
      index_type::level_type maxLevel;

      maxLevel = 0;

      for (size_t i = 0; i < list.size(); i++) {
        for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
          (list[i])->get(current_dim, curLevel, curIndex);

          if (curLevel > maxLevel) {
            maxLevel = curLevel;
          }
        }
      }

      return static_cast<size_t>(maxLevel);
    }

    void
    HashGridStorage::getLevelIndexMaskArraysForModEval(DataMatrix& level, DataMatrix& index,
        DataMatrix& mask, DataMatrix& offset) {
      index_type::level_type curLevel;
      index_type::level_type curIndex;

      for (size_t i = 0; i < list.size(); i++) {
        for (size_t current_dim = 0; current_dim < DIM; current_dim++) {
          (list[i])->get(current_dim, curLevel, curIndex);

          if (curLevel == 1) {
            level.set(i, current_dim, 0.0);
            index.set(i, current_dim, 0.0);
#if USE_DOUBLE_PRECISION
            uint64_t intmask = 0x0000000000000000;
#else
            uint32_t intmask = 0x00000000;
#endif
            mask.set(i, current_dim, *reinterpret_cast<float_t*>(&intmask));
            offset.set(i, current_dim, 1.0);
          } else if (curIndex == 1) {
            level.set(i, current_dim,
                      (-1.0) * static_cast<float_t>(1 << curLevel));
            index.set(i, current_dim, 0.0);
#if USE_DOUBLE_PRECISION
            uint64_t intmask = 0x0000000000000000;
#else
            uint32_t intmask = 0x00000000;
#endif
            mask.set(i, current_dim, *reinterpret_cast<float_t*>(&intmask));
            offset.set(i, current_dim, 2.0);
          } else if (curIndex
                     == static_cast<index_type::level_type>(((1 << curLevel) - 1))) {
            level.set(i, current_dim, static_cast<float_t>(1 << curLevel));
            index.set(i, current_dim, static_cast<float_t>(curIndex));
#if USE_DOUBLE_PRECISION
            uint64_t intmask = 0x0000000000000000;
#else
            uint32_t intmask = 0x00000000;
#endif
            mask.set(i, current_dim, *reinterpret_cast<float_t*>(&intmask));
            offset.set(i, current_dim, 1.0);
          } else {
            level.set(i, current_dim, static_cast<float_t>(1 << curLevel));
            index.set(i, current_dim, static_cast<float_t>(curIndex));
#if USE_DOUBLE_PRECISION
            uint64_t intmask = 0x8000000000000000;
#else
            uint32_t intmask = 0x80000000;
#endif
            mask.set(i, current_dim, *reinterpret_cast<float_t*>(&intmask));
            offset.set(i, current_dim, 1.0);
          }
        }
      }
    }

    void
    HashGridStorage::parseGridDescription(std::istream& istream) {
      int version;
      istream >> version;

      istream >> DIM;

      size_t num;
      istream >> num;

      // check whether grid was created with a version that is too new
      if (version > SERIALIZATION_VERSION) {
        if (version != 4) {
          std::ostringstream errstream;
          errstream << "Version of serialized grid (" << version
                    << ") is too new. Max. recognized version is "
                    << SERIALIZATION_VERSION << ".";
          throw generation_exception(errstream.str().c_str());
        }
      }

      //no bounding box, generate a trivial one
      if (version == 1 || version == 2) {
        // create a standard bounding box
        boundingBox = new BoundingBox(DIM);
      }

      // read the bounding box
      else if (version == 3 || version == 4) {
        // create a standard bounding box
        boundingBox = new BoundingBox(DIM);
        stretching = NULL;
        bUseStretching = false;
        DimensionBoundary tempBound;

        // reads the bounding box
        for (size_t i = 0; i < DIM; i++) {
          istream >> tempBound.leftBoundary;
          istream >> tempBound.rightBoundary;
          istream >> tempBound.bDirichletLeft;
          istream >> tempBound.bDirichletRight;

          boundingBox->setBoundary(i, tempBound);
        }
      } else if (version >= 5) {
        //      std::cout<<"Version 5 parse starts\n";
        int useStretching;
        DimensionBoundary tempBound;
        istream >> useStretching;

        if (useStretching == 0) {
          //BoundingBox

          // create a standard bounding box
          boundingBox = new BoundingBox(DIM);
          stretching = NULL;
          bUseStretching = false;

          // reads the boundary data
          for (size_t i = 0; i < DIM; i++) {
            istream >> tempBound.leftBoundary;
            istream >> tempBound.rightBoundary;
            istream >> tempBound.bDirichletLeft;
            istream >> tempBound.bDirichletRight;

            boundingBox->setBoundary(i, tempBound);
          }
        } else if (useStretching == 1) {
          //Stretching with analytic mode
          boundingBox = NULL;
          bUseStretching = true;
          Stretching1D* str1ds = new Stretching1D[DIM];
          DimensionBoundary* tempBounds = new DimensionBoundary[DIM];

          // reads the boundary data
          for (size_t i = 0; i < DIM; i++) {
            istream >> tempBounds[i].leftBoundary;
            istream >> tempBounds[i].rightBoundary;
            istream >> tempBounds[i].bDirichletLeft;
            istream >> tempBounds[i].bDirichletRight;
          }

          int stretchingType = 0;

          //Reads the 1D stretching data
          for (size_t i = 0; i < DIM; i++) {
            istream >> stretchingType;

            switch (stretchingType) {
              case 1:
                str1ds[i].type.assign("id");
                break;

              case 2:
                str1ds[i].type.assign("log");
                break;

              case 3:
                str1ds[i].type.assign("sinh");
                break;

              default:
                std::cout << "Stretching Type Unknown in parseGridDescription\n";
                break;
            }

            istream >> str1ds[i].x_0;
            istream >> str1ds[i].xsi;
          }

          stretching = new Stretching(DIM, tempBounds, str1ds);
          delete[] tempBounds;
          delete[] str1ds;
        } else if (useStretching == 2) {
          //Stretching with discrete Mode

          boundingBox = NULL;
          bUseStretching = true;

          // reads the boundary data, won't be used.
          for (size_t i = 0; i < DIM; i++) {
            istream >> tempBound.leftBoundary;
            istream >> tempBound.rightBoundary;
            istream >> tempBound.bDirichletLeft;
            istream >> tempBound.bDirichletRight;
          }

          int discreteLevel = 0;
          int vectorLength = 0;
          std::vector<float_t>* vec = new std::vector<float_t>[DIM];

          for (size_t i = 0; i < DIM; i++) {
            istream >> discreteLevel;
            vectorLength = static_cast<int>(pow(2.0, discreteLevel)) + 1;
            vec[i] = std::vector<float_t>(vectorLength, 0);

            for (int j = 0; j < vectorLength; j++) {
              istream >> vec[i][j];
            }
          }

          stretching = new Stretching(DIM, vec);
          delete[] vec;
        } else {
          std::cout << "Unknown Container Id Given in parseGridDescription\n";
        }
      }

      for (size_t i = 0; i < num; i++) {
        index_pointer index = new HashGridIndex(istream, version);
        list.push_back(index);
        map[index] = i;
      }

      // set's the grid point's leaf information which is not saved in version 1
      if (version == 1 || version == 4) {
        recalcLeafProperty();
      }
    }




  } // namespace base
} // namespace SGPP

