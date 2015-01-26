%rename(__getitem__) FullGrid::get(int index) const;
%rename(__setitem__) FullGrid::set(int index,double val);
%rename(__len__) FullGrid::getSize;
%rename(assign) FullGrid::operator=;
%include "std_vector.i"

%newobject SGPP::combigrid::FullGrid::createLinearFullGrid(size_t dim, vector<level_t> *inlevel);
%newobject SGPP::combigrid::FullGrid::createLinearBoundaryFullGrid(size_t dim, vector<level_t> *inlevel);


namespace std{
	%template(leveltvector) vector<SGPP::combigrid::FullGrid::level_t>;
}

using namespace SGPP::base;

namespace SGPP {
namespace combigrid{
class FullGrid{
public:
	  typedef GridStorage::index_type index_type;
	  typedef index_type::index_type index_t;
 	  typedef index_type::level_type level_t;
    FullGrid(const FullGrid &g);
    FullGrid(size_t indim, vector<level_t> *inlevel);
    FullGrid(size_t indim, size_t *inlevel);
    ~FullGrid();
    FullGrid operator=(const FullGrid &fg);
    size_t getSize();
    unsigned int *getLevel();
    unsigned int getLevel(size_t index);
    size_t getDim();     
    size_t length(size_t dim);
    double get(size_t index);
    void set(size_t index, double val);
    virtual double& at(size_t *index);
    virtual double getCoord(size_t d,size_t index);
    void fill(DataVector &v);  
    virtual void getCoords(size_t index,DataVector &v);
    std::string getCoordsString(size_t index);
    const char* getType();
    double eval(DataVector &v);
    double eval_for_frequency_domains(DataVector &v);
    BoundingBox* getBoundingBox();
    void setBoundingBox(BoundingBox *bBox);
/*static methods*/
   static FullGrid* createLinearFullGrid(size_t dim, vector<level_t> *inlevel);
   static FullGrid* createLinearBoundaryFullGrid(size_t dim, vector<level_t> *inlevel);


};
}
}

