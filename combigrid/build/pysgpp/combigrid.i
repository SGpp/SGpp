// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef SG_COMBIGRID

%include <std_shared_ptr.i>

%shared_ptr(sgpp::combigrid::CombigridOperation)
%shared_ptr(sgpp::combigrid::CombigridMultiOperation)
%shared_ptr(sgpp::combigrid::AbstractGrowthStrategy)
%shared_ptr(sgpp::combigrid::LinearGrowthStrategy)
%shared_ptr(sgpp::combigrid::ExponentialGrowthStrategy)
%shared_ptr(sgpp::combigrid::AbstractPointDistribution)
%shared_ptr(sgpp::combigrid::ClenshawCurtisDistribution)
%shared_ptr(sgpp::combigrid::UniformPointDistribution)
%shared_ptr(sgpp::combigrid::LejaPointDistribution)
%shared_ptr(sgpp::combigrid::AbstractPointOrdering)
%shared_ptr(sgpp::combigrid::ExponentialLevelorderPointOrdering)
%shared_ptr(sgpp::combigrid::IdentityPointOrdering)
%shared_ptr(sgpp::combigrid::AbstractPointHierarchy)
%shared_ptr(sgpp::combigrid::NestedPointHierarchy)
%shared_ptr(sgpp::combigrid::NonNestedPointHierarchy)
%shared_ptr(sgpp::combigrid::AbstractEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::BarycentricInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::BarycentricInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::LinearInterpolationEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::LinearInterpolationEvaluator>)
%shared_ptr(sgpp::combigrid::QuadratureEvaluator)
%shared_ptr(sgpp::combigrid::ArrayEvaluator<sgpp::combigrid::QuadratureEvaluator>)

%shared_ptr(sgpp::combigrid::AbstractCombigridStorage)
%shared_ptr(sgpp::combigrid::CombigridTreeStorage)

%shared_ptr(sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatArrayVector>)
%shared_ptr(sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatArrayVector>)

%shared_ptr(sgpp::combigrid::AbstractLevelEvaluator)
%shared_ptr(sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatScalarVector>)
%shared_ptr(sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatArrayVector>)


// %shared_ptr(sgpp::combigrid::AbstractLinearEvaluator<FloatScalarVector>)
// %shared_ptr(sgpp::combigrid::AbstractPermutationIterator)
// %shared_ptr(sgpp::combigrid::AbstractMultiStorage)
// %shared_ptr(sgpp::combigrid::AbstractMultiStorageIterator)



%include "combigrid/src/sgpp/combigrid/MultiFunction.hpp"
%include "combigrid/src/sgpp/combigrid/SingleFunction.hpp"
%include "combigrid/src/sgpp/combigrid/definitions.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/AbstractEvaluator.hpp"
%ignore sgpp::combigrid::FloatScalarVector::operator=;
%ignore sgpp::combigrid::FloatArrayVector::operator=;
%ignore sgpp::combigrid::FloatArrayVector::operator[];
%include "combigrid/src/sgpp/combigrid/algebraic/FloatScalarVector.hpp"
%include "combigrid/src/sgpp/combigrid/algebraic/FloatArrayVector.hpp"
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

namespace sgpp {
namespace combigrid {
    %template(AbstractEvaluator_FloatScalarVector) AbstractEvaluator<FloatScalarVector>;
    %template(AbstractEvaluator_FloatArrayVector) AbstractEvaluator<FloatArrayVector>;
    %template(AbstractLinearEvaluator_FloatScalarVector) AbstractLinearEvaluator<FloatScalarVector>;
    %template(AbstractLinearEvaluator_FloatArrayVector) AbstractLinearEvaluator<FloatArrayVector>;
    %template(FloatArrayVectorMultiStorage) AbstractMultiStorage<FloatArrayVector>;
    %template(FloatScalarVectorMultiStorage) AbstractMultiStorage<FloatScalarVector>;
    %template(FloatArrayVectorMultiStorageIterator) AbstractMultiStorageIterator<FloatArrayVector>;
    %template(FloatScalarVectorMultiStorageIterator) AbstractMultiStorageIterator<FloatScalarVector>;

    // %template(AbstractMultiStorage_uint8_t) AbstractMultiStorage<uint8_t>;
    // %template(TreeStorage_uint8_t) TreeStorage<uint8_t>;
}
}

%include "combigrid/src/sgpp/combigrid/operation/multidim/FullGridTensorEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/AdaptiveRefinementStrategy.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/AbstractLevelEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/BarycentricInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp"
%include "combigrid/src/sgpp/combigrid/operation/onedim/ArrayEvaluator.hpp"

// %include "combigrid/src/sgpp/combigrid/serialization/AbstractSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/DefaultSerializationStrategy.hpp"
// %include "combigrid/src/sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp"

namespace sgpp{
namespace combigrid {
    %template(ScalarAbstractFullGridEvaluator) AbstractFullGridEvaluator<FloatScalarVector>;
    %template(ArrayAbstractFullGridEvaluator) sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatArrayVector>;

    %template(ScalarFullGridTensorEvaluator) FullGridTensorEvaluator<FloatScalarVector>;
    %template(ArrayFullGridTensorEvaluator) sgpp::combigrid::FullGridTensorEvaluator<sgpp::combigrid::FloatArrayVector>;
    %template(ScalarCombigridEvaluator) CombigridEvaluator<FloatScalarVector>;
    %template(ArrayCombigridEvaluator) CombigridEvaluator<FloatArrayVector>;

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

    %template(FloatScalarAbstractLinearEvaluatorVector) vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatScalarVector>>>;
    %template(FloatArrayAbstractLinearEvaluatorVector) vector<std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>;
    %template(AbstractPointHierarchyVector) vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>>;

    %template(FloatScalarVectorVector) vector<sgpp::combigrid::FloatScalarVector>;
    %template(FloatArrayVectorVector) vector<sgpp::combigrid::FloatArrayVector>;
    // %template(MultidimFunction) std::function<double(sgpp::base::DataVector const &)>;
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
namespace sgpp
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
static double PythonCallBackFunc(int len, double* a, void *clientdata)
{
  PyObject *func, *lst, *arglist;
   PyObject *result;
   double    dres = 0;

   // get Python function
   func = (PyObject *) clientdata;
   // build argument list (only Python list, convert double* to Python list)
   lst = PyList_New(len);        // alternatively: PyList_New(len)
   if (!lst) {
     PyErr_SetString(PyExc_TypeError, "No data provided!");
     return NULL;
   }
   for (int i=0; i<len; i++) {
     // create new Python double
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

static sgpp::combigrid::MultiFunction multiFunc(PyObject *func) {
    return MultiFunction(PyFuncWrapper<double, base::DataVector const &>(func, [=](base::DataVector const &vec) -> double {
        base::DataVector v = vec;
        return PythonCallBackFunc(vec.getSize(), &v[0], (void *)func);
    }));
}

/* This function has to match the type and parameters of the
 * C++ callback functions that are used. However, the
 * clientdata pointer is used for holding a reference to a
 * Python callable object.
 */
static double PythonSingleCallBackFunc(double a, void *clientdata)
{
  PyObject *func, *arglist;
   PyObject *result;
   double    dres = 0;

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

static sgpp::combigrid::SingleFunction singleFunc(PyObject *func) {
    return SingleFunction(PyFuncWrapper<double, double>(func, [=](double x) -> double {
        return PythonSingleCallBackFunc(x, (void *)func);
    }));
}

}
}
%}

#endif
