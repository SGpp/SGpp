// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef SG_COMBIGRID

%include <std_shared_ptr.i>

%shared_ptr(SGPP::combigrid::CombigridOperation)
%shared_ptr(SGPP::combigrid::CombigridMultiOperation)
%shared_ptr(SGPP::combigrid::AbstractGrowthStrategy)
%shared_ptr(SGPP::combigrid::LinearGrowthStrategy)
%shared_ptr(SGPP::combigrid::ExponentialGrowthStrategy)
%shared_ptr(SGPP::combigrid::AbstractPointDistribution)
%shared_ptr(SGPP::combigrid::ClenshawCurtisDistribution)
%shared_ptr(SGPP::combigrid::UniformPointDistribution)
%shared_ptr(SGPP::combigrid::LejaPointDistribution)
%shared_ptr(SGPP::combigrid::AbstractPointOrdering)
%shared_ptr(SGPP::combigrid::ExponentialLevelorderPointOrdering)
%shared_ptr(SGPP::combigrid::IdentityPointOrdering)
%shared_ptr(SGPP::combigrid::AbstractPointHierarchy)
%shared_ptr(SGPP::combigrid::NestedPointHierarchy)
%shared_ptr(SGPP::combigrid::NonNestedPointHierarchy)
%shared_ptr(SGPP::combigrid::AbstractEvaluator<SGPP::combigrid::ScalarVector<SGPP::float_t>>)
%shared_ptr(SGPP::combigrid::AbstractEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>)
%shared_ptr(SGPP::combigrid::AbstractLinearEvaluator<SGPP::combigrid::ScalarVector<SGPP::float_t>>)
%shared_ptr(SGPP::combigrid::AbstractLinearEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>)
%shared_ptr(SGPP::combigrid::BarycentricInterpolationEvaluator)
%shared_ptr(SGPP::combigrid::ArrayEvaluator<SGPP::combigrid::BarycentricInterpolationEvaluator>)
%shared_ptr(SGPP::combigrid::LinearInterpolationEvaluator)
%shared_ptr(SGPP::combigrid::ArrayEvaluator<SGPP::combigrid::LinearInterpolationEvaluator>)
%shared_ptr(SGPP::combigrid::QuadratureEvaluator)
%shared_ptr(SGPP::combigrid::ArrayEvaluator<SGPP::combigrid::QuadratureEvaluator>)

%shared_ptr(SGPP::combigrid::AbstractCombigridStorage)
%shared_ptr(SGPP::combigrid::CombigridTreeStorage)

%shared_ptr(SGPP::combigrid::AbstractFullGridEvaluator<SGPP::combigrid::ScalarVector<SGPP::float_t>>)
%shared_ptr(SGPP::combigrid::AbstractFullGridEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>)
%shared_ptr(SGPP::combigrid::FullGridTensorEvaluator<SGPP::combigrid::ScalarVector<SGPP::float_t>>)
%shared_ptr(SGPP::combigrid::FullGridTensorEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>)

%shared_ptr(SGPP::combigrid::CombigridEvaluator<SGPP::combigrid::ScalarVector<SGPP::float_t>>)
%shared_ptr(SGPP::combigrid::CombigridEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>)


// %shared_ptr(SGPP::combigrid::AbstractLinearEvaluator<ScalarVector<SGPP::float_t>>)
// %shared_ptr(SGPP::combigrid::AbstractPermutationIterator)
// %shared_ptr(SGPP::combigrid::AbstractMultiStorage)
// %shared_ptr(SGPP::combigrid::AbstractMultiStorageIterator)



%include "combigrid/src/sgpp/combigrid/MultiFunction.hpp"
%include "combigrid/src/sgpp/combigrid/SingleFunction.hpp"
%include "combigrid/src/sgpp/combigrid/definitions.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractEvaluator.hpp"
%ignore SGPP::combigrid::ScalarVector::operator=;
%ignore SGPP::combigrid::ArrayVector::operator=;
%ignore SGPP::combigrid::ArrayVector::operator[];
%include "combigrid/src/sgpp/combigrid/algebraic/ScalarVector.hpp"
%include "combigrid/src/sgpp/combigrid/algebraic/ArrayVector.hpp"
%include "combigrid/src/sgpp/combigrid/common/MultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/AbstractPermutationIterator.hpp"
%include "combigrid/src/sgpp/combigrid/storage/IterationPolicy.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorage.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractCombigridStorage.hpp"
%include "combigrid/src/sgpp/combigrid/storage/FunctionLookupTable.hpp"
%include "combigrid/src/sgpp/combigrid/storage/tree/TreeStorage.hpp"

%include "combigrid/src/sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/ClenshawCurtisDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/ExponentialGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/ExponentialLevelorderPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp"

%include "combigrid/src/sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp"

%include "combigrid/src/sgpp/combigrid/operation/multidim/AbstractFullGridEvaluator.hpp"

%include "combigrid/src/sgpp/combigrid/threading/ThreadPool.hpp"

namespace SGPP {
namespace combigrid {
    %template(PyFloatScalarVector) ScalarVector<float_t>;
    %template(PyFloatArrayVector) ArrayVector<float_t, ScalarVector<float_t>>;
    %template(AbstractEvaluator_FloatScalarVector) AbstractEvaluator<ScalarVector<float_t>>;
    %template(AbstractEvaluator_FloatArrayVector) AbstractEvaluator<ArrayVector<float_t, ScalarVector<float_t>>>;
    %template(AbstractLinearEvaluator_FloatScalarVector) AbstractLinearEvaluator<ScalarVector<float_t>>;
    %template(AbstractLinearEvaluator_FloatArrayVector) AbstractLinearEvaluator<ArrayVector<float_t, ScalarVector<float_t>>>;
    %template(FloatArrayVectorMultiStorage) AbstractMultiStorage<ArrayVector<float_t, ScalarVector<float_t>>>;
    %template(FloatScalarVectorMultiStorage) AbstractMultiStorage<ScalarVector<float_t>>;
    %template(FloatArrayVectorMultiStorageIterator) AbstractMultiStorageIterator<ArrayVector<float_t, ScalarVector<float_t>>>;
    %template(FloatScalarVectorMultiStorageIterator) AbstractMultiStorageIterator<ScalarVector<float_t>>;
    
    // %template(AbstractMultiStorage_uint8_t) AbstractMultiStorage<uint8_t>;
    // %template(TreeStorage_uint8_t) TreeStorage<uint8_t>;
}
}

%include "combigrid/src/sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp"

// %include "combigrid/src/sgpp/combigrid/serialization/AbstractSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/DefaultSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp"

namespace SGPP {
namespace combigrid {
    %template(ScalarAbstractFullGridEvaluator) AbstractFullGridEvaluator<ScalarVector<float_t>>;
    %template(ArrayAbstractFullGridEvaluator) SGPP::combigrid::AbstractFullGridEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>;

    %template(ScalarFullGridTensorEvaluator) FullGridTensorEvaluator<ScalarVector<float_t>>;
    %template(ArrayFullGridTensorEvaluator) SGPP::combigrid::FullGridTensorEvaluator<SGPP::combigrid::ArrayVector<SGPP::float_t, SGPP::combigrid::ScalarVector<SGPP::float_t>>>;
    %template(ScalarCombigridEvaluator) CombigridEvaluator<ScalarVector<float_t>>;
    %template(ArrayCombigridEvaluator) CombigridEvaluator<ArrayVector<float_t, ScalarVector<float_t>>>;
    
    %template(ArrayBarycentricInterpolationEvaluator) ArrayEvaluator<BarycentricInterpolationEvaluator>;
    %template(ArrayLinearInterpolationEvaluator) ArrayEvaluator<LinearInterpolationEvaluator>;
    %template(ArrayQuadratureEvaluator) ArrayEvaluator<QuadratureEvaluator>;
    
    // %template(AbstractSerializationStrategy_uint8_t) AbstractSerializationStrategy<std::shared_ptr<TreeStorage<uint8_t>>>;
    // %template(AbstractSerializationStrategy_uint8_t) AbstractSerializationStrategy<uint8_t>;
    // %template(DefaultSerializationStrategy_uint8_t) DefaultSerializationStrategy<uint8_t>;
    // %template(LevelStructureSerializationStrategy) TreeStorageSerializationStrategy<uint8_t>;
    
}
}

namespace std {

    %template(FloatScalarAbstractLinearEvaluatorVector) vector<std::shared_ptr<SGPP::combigrid::AbstractLinearEvaluator<SGPP::combigrid::ScalarVector<float_t>>>>;
    %template(FloatArrayAbstractLinearEvaluatorVector) vector<std::shared_ptr<SGPP::combigrid::AbstractLinearEvaluator<SGPP::combigrid::ArrayVector<float_t, SGPP::combigrid::ScalarVector<float_t>>>>>;
    %template(AbstractPointHierarchyVector) vector<std::shared_ptr<SGPP::combigrid::AbstractPointHierarchy>>;
    
    %template(FloatScalarVectorVector) vector<SGPP::combigrid::FloatScalarVector>;
    %template(FloatArrayVectorVector) vector<SGPP::combigrid::FloatArrayVector>;
    // %template(MultidimFunction) std::function<SGPP::float_t(SGPP::base::DataVector const &)>;
}

%include "combigrid/src/sgpp/combigrid/common/AbstractPermutationIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/MultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/common/BoundedSumMultiIndexIterator.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp"
%include "combigrid/src/sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp"
%include "combigrid/src/sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp"
%include "combigrid/src/sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp"

%include "combigrid/src/sgpp/combigrid/numeric/KahanAdder.hpp"
%include "combigrid/src/sgpp/combigrid/storage/AbstractCombigridStorage.hpp"
%include "combigrid/src/sgpp/combigrid/operation/CombigridOperation.hpp"
%include "combigrid/src/sgpp/combigrid/operation/CombigridMultiOperation.hpp"

%inline %{
namespace SGPP
{
namespace combigrid
{

#include <functional>

template<typename Out, typename In> class PyFuncWrapper {
    PyObject *func;
    std::function<Out(In)> stdfunc;
    
public:
    PyFuncWrapper(PyObject *func, std::function<Out(In)> stdfunc) : func(func), stdfunc(stdfunc) {
        Py_INCREF(func);
    }
    
    PyFuncWrapper(PyFuncWrapper<Out, In> const &other) : func(other.func), stdfunc(other.stdfunc) {
        Py_INCREF(func);
    }
    
    ~PyFuncWrapper() {
        Py_DECREF(func);
    }
    
    Out operator()(In param) {
        return stdfunc(param);
    }
};

/* This function has to match the type and parameters of the
 * C++ callback functions that are used. However, the
 * clientdata pointer is used for holding a reference to a 
 * Python callable object. 
 */
static float_t PythonCallBackFunc(int len, float_t* a, void *clientdata)
{
  PyObject *func, *lst, *arglist;
   PyObject *result;
   float_t    dres = 0;

   // get Python function
   func = (PyObject *) clientdata;
   // build argument list (only Python list, convert float_t* to Python list)
   lst = PyList_New(len);        // alternatively: PyList_New(len)
   if (!lst) {
     PyErr_SetString(PyExc_TypeError, "No data provided!");
     return NULL;
   }
   for (int i=0; i<len; i++) {
     // create new Python float_t
     PyObject *num = PyFloat_FromDouble(a[i]);
     if (!num) {
       PyErr_SetString(PyExc_TypeError, "No data in list!");
       Py_DECREF(lst);
       return NULL;
     }
     // steals reference to num:
     PyList_SetItem(lst, i, num); // alternatively: PyList_SET_ITEM()
   }
   // build list of one Python object
   arglist = Py_BuildValue("(O)", lst);
   // call Python
   result = PyEval_CallObject(func,arglist);
   // trash arglist and lst
   Py_DECREF(arglist);
   Py_DECREF(lst);
   if (result) {
     dres = PyFloat_AsDouble(result);
   }
   Py_XDECREF(result);
   
   return dres;
}

static SGPP::combigrid::MultiFunction multiFunc(PyObject *func) {
    return MultiFunction(PyFuncWrapper<SGPP::float_t, base::DataVector const &>(func, [=](base::DataVector const &vec) -> SGPP::float_t {
        base::DataVector v = vec;
        return PythonCallBackFunc(vec.getSize(), &v[0], (void *)func);
    }));
}

/* This function has to match the type and parameters of the
 * C++ callback functions that are used. However, the
 * clientdata pointer is used for holding a reference to a 
 * Python callable object. 
 */
static float_t PythonSingleCallBackFunc(float_t a, void *clientdata)
{
  PyObject *func, *arglist;
   PyObject *result;
   float_t    dres = 0;
   
   // std::cout << "Call single func with parameter " << a << "\n";

   // get Python function
   func = (PyObject *) clientdata;
   
   PyObject *param = PyFloat_FromDouble(a);
   if (!param) {
        PyErr_SetString(PyExc_TypeError, "Conversion error in singleFunc()");
   }
   
   // build list of one Python object
   arglist = Py_BuildValue("(O)", param);
   // call Python
   result = PyEval_CallObject(func,arglist);
   // trash arglist
   Py_DECREF(arglist);
   Py_DECREF(param);
   if (result) {
     dres = PyFloat_AsDouble(result);
   }
   Py_XDECREF(result);
   // std::cout << "Result: " << dres << "\n";
   return dres;
}

static SGPP::combigrid::SingleFunction singleFunc(PyObject *func) {
    return SingleFunction(PyFuncWrapper<SGPP::float_t, SGPP::float_t>(func, [=](SGPP::float_t x) -> SGPP::float_t {
        return PythonSingleCallBackFunc(x, (void *)func);
    }));
}

}
}
%}

#endif