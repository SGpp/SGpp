/*
 * SimpleMultiClassGenerator.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: katrin
 */

#ifndef SIMPLEMULTICLASSGENERATOR_HPP
#define SIMPLEMULTICLASSGENERATOR_HPP


#include <sgpp/base/grid/Grid.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
class SimpleMultiClassGenerator {
public:
    SimpleMultiClassGenerator(int dim, int classes, int level);
    virtual ~SimpleMultiClassGenerator();
    
    double getEval(int classID, int pointAt);
    
    std::vector<base::Grid*> getMulitGrid();
    
private:
    std::vector<base::Grid*> grids;
};


} /* namespace datadriven */
} /* namespace sgpp */
#endif /* SIMPLEMULTICLASSGENERATOR_HPP */
