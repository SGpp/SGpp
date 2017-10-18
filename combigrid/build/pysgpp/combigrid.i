// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef SG_COMBIGRID

%include <std_shared_ptr.i>

%shared_ptr(sgpp::combigrid::CombigridOperation)
%shared_ptr(sgpp::combigrid::CombigridMultiOperation)
%shared_ptr(sgpp::combigrid::CombigridTensorOperation)
%shared_ptr(sgpp::combigrid::AbstractGrowthStrategy)
%shared_ptr(sgpp::combigrid::LinearGrowthStrategy)
%shared_ptr(sgpp::combigrid::ExponentialGrowthStrategy)
%shared_ptr(sgpp::combigrid::AbstractPointDistribution)
%shared_ptr(sgpp::combigrid::ClenshawCurtisDistribution)
%shared_ptr(sgpp::combigrid::UniformPointDistribution)
%shared_ptr(sgpp::combigrid::UniformBoundaryPointDistribution)
%shared_ptr(sgpp::combigrid::ChebyshevDistribution)
%shared_ptr(sgpp::combigrid::LejaPointDistribution)
%shared_ptr(sgpp::combigrid::L2LejaPointDistribution)
%shared_ptr(sgpp::combigrid::AbstractPointOrdering)
%shared_ptr(sgpp::combigrid::ExponentialLevelorderPointOrdering)
%shared_ptr(sgpp::combigrid::IdentityPointOrdering)
%shared_ptr(sgpp::combigrid::AbstractPointHierarchy)
%shared_ptr(sgpp::combigrid::NestedPointHierarchy)
%shared_ptr(sgpp::combigrid::NonNestedPointHierarchy)
%shared_ptr(sgpp::combigrid::AbstractEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::AbstractEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::PolynomialInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::PolynomialInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::LinearInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::LinearInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::CubicSplineInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::CubicSplineInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::QuadratureEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::QuadratureEvaluator>)
%shared_ptr(sgpp::combigrid::BSplineInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::BSplineInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::InterpolationCoefficientEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::InterpolationCoefficientEvaluator>)

%shared_ptr(sgpp::combigrid::AbstractCombigridStorage)
%shared_ptr(sgpp::combigrid::CombigridTreeStorage)
%shared_ptr(sgpp::combigrid::AbstractMultiStorage<double>)
%shared_ptr(sgpp::combigrid::AbstractMultiStorage<uint8_t>)
%shared_ptr(sgpp::combigrid::TreeStorage<double>)
%shared_ptr(sgpp::combigrid::TreeStorage<uint8_t>)

%shared_ptr(sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridLinearEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridLinearEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridLinearEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearCallbackEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearCallbackEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearCallbackEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearGridBasedEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearGridBasedEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::FullGridLinearGridBasedEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatTensorVector>)

%shared_ptr(sgpp::combigrid::AbstractLevelEvaluator)
%shared_ptr(sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatTensorVector>)
%shared_ptr(sgpp::combigrid::CombigridTensorEvaluator<sgpp::combigrid::FloatTensorVector>)


%shared_ptr(sgpp::combigrid::AbstractInfiniteFunctionBasis1D)
%shared_ptr(sgpp::combigrid::OrthogonalPolynomialBasis1D)
%shared_ptr(sgpp::combigrid::MonomialFunctionBasis1D)

%shared_ptr(sgpp::combigrid::LevelManager)
%shared_ptr(sgpp::combigrid::AveragingLevelManager)
%shared_ptr(sgpp::combigrid::WeightedRatioLevelManager)
%shared_ptr(sgpp::combigrid::VarianceLevelManager)
%shared_ptr(sgpp::combigrid::RegularLevelManager)

%shared_ptr(sgpp::combigrid::TensorGrid)
%shared_ptr(sgpp::combigrid::ThreadPool)

%shared_ptr(std::mutex)


// %shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<FloatScalarVector>)
// %shared_ptr(sgpp::combigrid::AbstractPermutationIterator)
// %shared_ptr(sgpp::combigrid::AbstractMultiStorage)
// %shared_ptr(sgpp::combigrid::AbstractMultiStorageIterator)



%include "combigrid/src/sgpp/combigrid/GeneralFunction.hpp"
%include "combigrid/src/sgpp/combigrid/threading/ThreadPool.hpp"

namespace sgpp {
namespace combigrid {

    %template(PyMultiFunction) GeneralFunction<double, base::DataVector const &>;
    %template(PySingleFunction) GeneralFunction<double, double>;
    %template(PyTask) GeneralFunction1<void>;
    %template(PyIdleFunction) GeneralFunction<void, ThreadPool &>;
}
}


%include "combigrid/src/sgpp/combigrid/definitions.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractEvaluator.hpp"

%include "combigrid/src/sgpp/combigrid/storage/tree/TreeStorage.hpp"

%ignore sgpp::combigrid::FloatScalarVector::operator=;
%ignore sgpp::combigrid::FloatArrayVector::operator=;
%ignore sgpp::combigrid::FloatArrayVector::operator[];
%ignore sgpp::combigrid::FloatTensorVector::operator=;
%ignore sgpp::combigrid::FloatTensorVector::operator[];

%ignore sgpp::combigrid::operator<<;
%ignore sgpp::combigrid::operator>>;

%include "combigrid/src/sgpp/combigrid/algebraic/FloatScalarVector.hpp"
%include "combigrid/src/sgpp/combigrid/algebraic/FloatArrayVector.hpp"
%include "combigrid/src/sgpp/combigrid/algebraic/FloatTensorVector.hpp"
%include "combigrid/src/sgpp/combigrid/common/MultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/AbstractPermutationIterator.hpp"
%include "combigrid/src/sgpp/combigrid/storage/IterationPolicy.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorage.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractCombigridStorage.hpp"
%include "combigrid/src/sgpp/combigrid/storage/FunctionLookupTable.hpp"

%include "combigrid/src/sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp"
%include "combigrid/src/sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp"
%include "combigrid/src/sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp"

%include "combigrid/src/sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/L2LejaPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/UniformBoundaryPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/ChebyshevDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/TensorGrid.hpp"

%include "combigrid/src/sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp"

%include "combigrid/src/sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridLinearEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearCallbackEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearGridBasedEvaluator.hpp"

%include "combigrid/src/sgpp/combigrid/threading/ThreadPool.hpp"

namespace sgpp {
namespace combigrid {
    %template(AbstractEvaluator_FloatScalarVector) AbstractEvaluator<FloatScalarVector>;
    %template(AbstractEvaluator_FloatArrayVector) AbstractEvaluator<FloatArrayVector>;
    %template(AbstractEvaluator_FloatTensorVector) AbstractEvaluator<FloatTensorVector>;
    %template(AbstractLinearEvaluator_FloatScalarVector) AbstractLinearEvaluator<FloatScalarVector>;
    %template(AbstractLinearEvaluator_FloatArrayVector) AbstractLinearEvaluator<FloatArrayVector>;
    %template(AbstractLinearEvaluator_FloatTensorVector) AbstractLinearEvaluator<FloatTensorVector>;
    %template(FloatScalarVectorMultiStorage) AbstractMultiStorage<FloatScalarVector>;
    %template(FloatArrayVectorMultiStorage) AbstractMultiStorage<FloatArrayVector>;
    %template(FloatTensorVectorMultiStorage) AbstractMultiStorage<FloatTensorVector>;
    %template(FloatScalarVectorMultiStorageIterator) AbstractMultiStorageIterator<FloatScalarVector>;
    %template(FloatArrayVectorMultiStorageIterator) AbstractMultiStorageIterator<FloatArrayVector>;
    %template(FloatTensorVectorMultiStorageIterator) AbstractMultiStorageIterator<FloatTensorVector>;

    %template(DoubleAbstractMultiStorage) AbstractMultiStorage<double>;
    %template(Uint8AbstractMultiStorage) AbstractMultiStorage<uint8_t>;
    %template(DoubleTreeStorage) TreeStorage<double>;
    %template(Uint8TreeStorage) TreeStorage<uint8_t>;
    %template(PyGridFunction) GeneralFunction<std::shared_ptr<TreeStorage<double>>, std::shared_ptr<TensorGrid>>;
}
}

%include "combigrid/src/sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/PolynomialInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/CubicSplineInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/InterpolationCoefficientEvaluator.hpp"

%include "combigrid/src/sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/Configurations.hpp"

// %include "combigrid/src/sgpp/combigrid/serialization/AbstractSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/DefaultSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp"

namespace sgpp{
namespace combigrid {
    %template(ScalarAbstractFullGridEvaluator) AbstractFullGridEvaluator<FloatScalarVector>;
    %template(ArrayAbstractFullGridEvaluator) AbstractFullGridEvaluator<FloatArrayVector>;
    %template(TensorAbstractFullGridEvaluator) AbstractFullGridEvaluator<FloatTensorVector>;

    %template(ScalarAbstractFullGridLinearEvaluator) AbstractFullGridLinearEvaluator<FloatScalarVector>;
    %template(ArrayAbstractFullGridLinearEvaluator) AbstractFullGridLinearEvaluator<FloatArrayVector>;
    %template(TensorAbstractFullGridLinearEvaluator) AbstractFullGridLinearEvaluator<FloatTensorVector>;

    %template(ScalarFullGridLinearCallbackEvaluator) FullGridLinearCallbackEvaluator<FloatScalarVector>;
    %template(ArrayFullGridLinearCallbackEvaluator) FullGridLinearCallbackEvaluator<FloatArrayVector>;
    %template(TensorFullGridLinearCallbackEvaluator) FullGridLinearCallbackEvaluator<FloatTensorVector>;

    %template(ScalarFullGridLinearGridBasedEvaluator) FullGridLinearGridBasedEvaluator<FloatScalarVector>;
    %template(ArrayFullGridLinearGridBasedEvaluator) FullGridLinearGridBasedEvaluator<FloatArrayVector>;
    %template(TensorFullGridLinearGridBasedEvaluator) FullGridLinearGridBasedEvaluator<FloatTensorVector>;

    %template(ScalarCombigridEvaluator) CombigridEvaluator<FloatScalarVector>;
    %template(ArrayCombigridEvaluator) CombigridEvaluator<FloatArrayVector>;
    %template(TensorCombigridEvaluator) CombigridEvaluator<FloatTensorVector>;

    %template(ArrayPolynomialInterpolationEvaluator) ArrayEvaluator<PolynomialInterpolationEvaluator>;
    %template(ArrayLinearInterpolationEvaluator) ArrayEvaluator<LinearInterpolationEvaluator>;
    %template(ArrayCubicSplineInterpolationEvaluator) ArrayEvaluator<CubicSplineInterpolationEvaluator>;
    %template(ArrayQuadratureEvaluator) ArrayEvaluator<QuadratureEvaluator>;
    %template(ArrayBSplineInterpolationEvaluator) ArrayEvaluator<BSplineInterpolationEvaluator>;

    // %template(AbstractSerializationStrategy_uint8_t) AbstractSerializationStrategy<std::shared_ptr<TreeStorage<std::uint8_t>>>;
    // %template(AbstractSerializationStrategy_uint8_t) AbstractSerializationStrategy<std::uint8_t>;
    // %template(DefaultSerializationStrategy_uint8_t) DefaultSerializationStrategy<std::uint8_t>;
    // %template(LevelStructureSerializationStrategy) TreeStorageSerializationStrategy<std::uint8_t>;

}
}

namespace std {

    %template(FloatScalarAbstractLinearEvaluatorVector) vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatScalarVector>>>;
    %template(FloatArrayAbstractLinearEvaluatorVector) vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>;
    %template(FloatTensorAbstractLinearEvaluatorVector) vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatTensorVector>>>;
    %template(AbstractPointHierarchyVector) vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>>;
    %template(OrthogonalPolynomialBasis1DVector) std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>;
    %template(OrthogonalPolynomialBasisTypeVector) std::vector<sgpp::combigrid::OrthogonalPolynomialBasisType>;

    %template(FloatScalarVectorVector) vector<sgpp::combigrid::FloatScalarVector>;
    %template(FloatArrayVectorVector) vector<sgpp::combigrid::FloatArrayVector>;
    %template(FloatTensorVectorVector) vector<sgpp::combigrid::FloatTensorVector>;
    %template(DataVectorVector) vector<sgpp::base::DataVector>;

    // %template(PyTaskVector) std::vector<sgpp::combigrid::GeneralFunction1<void>>;

    // %template(CombiHierarchiesCollection) std::vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>>;
    // %template(CombiEvaluatorsCollection) std::vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatScalarVector>>>;
    // %template(CombiEvaluatorsMultiCollection) std::vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>;
    // %template(MultidimFunction) std::function<double(sgpp::base::DataVector const &)>;
}

%include "combigrid/src/sgpp/combigrid/common/AbstractPermutationIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/MultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/GridConversion.hpp"

%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp"

%include "combigrid/src/sgpp/combigrid/numeric/KahanAdder.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractCombigridStorage.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/LevelManager.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/RegularLevelManager.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/VarianceLevelManager.hpp"
%include "combigrid/src/sgpp/combigrid/operation/CombigridOperation.hpp"
%include "combigrid/src/sgpp/combigrid/operation/CombigridMultiOperation.hpp"
%include "combigrid/src/sgpp/combigrid/operation/CombigridTensorOperation.hpp"

%include "combigrid/src/sgpp/combigrid/threading/ThreadPool.hpp"
%include "combigrid/src/sgpp/combigrid/threading/PtrGuard.hpp"

// %include "combigrid/src/sgpp/combigrid/utils/BinaryHeap.hpp" // is a template
%include "combigrid/src/sgpp/combigrid/utils/Stopwatch.hpp"
%include "combigrid/src/sgpp/combigrid/utils/Utils.hpp"

// experimental

%feature("director") sgpp::combigrid::GeneralFunctionDirector;
%feature("director") sgpp::combigrid::GeneralFunctionDirector1;
%include "combigrid/src/sgpp/combigrid/GeneralFunctionDirector.hpp"

namespace sgpp {
namespace combigrid {
    %template(MultiFunctionDirector) GeneralFunctionDirector<double, base::DataVector const &>;
    %template(SingleFunctionDirector) GeneralFunctionDirector<double, double>;
    %template(GridFunctionDirector) GeneralFunctionDirector<std::shared_ptr<TreeStorage<double>>, std::shared_ptr<TensorGrid>>;
    %template(ThreadPoolTaskDirector) GeneralFunctionDirector1<void>;
    %template(ThreadPoolIdleCallbackDirector) GeneralFunctionDirector<void, ThreadPool &>;
}
}


%pythoncode %{
class MFDirectorImpl(MultiFunctionDirector):
    def __init__(self):
        super(MFDirectorImpl, self).__init__()

    def setFuncObj(self, funcObj):
        self.funcObj = funcObj

    def eval(self, vec):
        return self.funcObj(vec)

def multiFunc(funcObj):
    dir = MFDirectorImpl()
    dir.setFuncObj(funcObj)
    mf = dir.toFunction()
    dir.__disown__()
    return mf

class SFDirectorImpl(SingleFunctionDirector):
    def __init__(self):
        super(SFDirectorImpl, self).__init__()

    def setFuncObj(self, funcObj):
        self.funcObj = funcObj

    def eval(self, vec):
        return self.funcObj(vec)

def singleFunc(funcObj):
    dir = SFDirectorImpl()
    dir.setFuncObj(funcObj)
    mf = dir.toFunction()
    dir.__disown__()
    return mf

class GFDirectorImpl(GridFunctionDirector):
    def __init__(self):
        super(GFDirectorImpl, self).__init__()

    def setFuncObj(self, funcObj):
        self.funcObj = funcObj

    def eval(self, vec):
        return self.funcObj(vec)

def gridFunc(funcObj):
    dir = GFDirectorImpl()
    dir.setFuncObj(funcObj)
    mf = dir.toFunction()
    dir.__disown__()
    return mf

class TaskDirectorImpl(ThreadPoolTaskDirector):
    def __init__(self):
        super(TaskDirectorImpl, self).__init__()

    def setFuncObj(self, funcObj):
        self.funcObj = funcObj

    def eval(self):
        return self.funcObj()

def taskFunc(funcObj):
    dir = TaskDirectorImpl()
    dir.setFuncObj(funcObj)
    f = dir.toFunction()
    dir.__disown__()
    return f

class IdleCallbackDirectorImpl(ThreadPoolIdleCallbackDirector):
    def __init__(self):
        super(IdleCallbackDirectorImpl, self).__init__()

    def setFuncObj(self, funcObj):
        self.funcObj = funcObj

    def eval(self):
        return self.funcObj()

def idleCallbackFunc(funcObj):
    dir = IdleCallbackDirectorImpl()
    dir.setFuncObj(funcObj)
    f = dir.toFunction()
    dir.__disown__()
    return f
%}

// does some exception handling according to the SWIG website
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

#endif
