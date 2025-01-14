// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// %include "combigrid/src/sgpp/combigrid/tools/IndexVectorRange.hpp"
// make IndexVectorRange iterable from Python
// https://stackoverflow.com/questions/35291053/how-to-make-a-c-class-iterable-from-python-using-swig
%inline %{
namespace sgpp{
namespace combigrid{
class StopIterator {};
}
}
%}

// remove SWIG warnings
%ignore sgpp::combigrid::IndexVectorIterator::operator=;
%ignore sgpp::combigrid::IndexVectorIterator::operator[];
%ignore sgpp::combigrid::IndexVectorIterator::operator++;
%ignore sgpp::combigrid::IndexVectorIterator::operator--;
%ignore sgpp::combigrid::IndexVectorIterator::operator+;
%ignore sgpp::combigrid::LevelVectorTools::Hash;

// get rid of warning 503
// %ignore sgpp::combigrid::operator+;

// get rid of warning 401
// why does none of this work?!?!
// %include <std/std_container.i>
// %ignore std::iterator< std::random_access_iterator_tag,sgpp::combigrid::IndexVector,size_t,sgpp::combigrid::IndexVector *,sgpp::combigrid::IndexVector & >;
// %template(IndexVectorIterator_base) std::iterator< std::random_access_iterator_tag,sgpp::combigrid::IndexVector,size_t,sgpp::combigrid::IndexVector *,sgpp::combigrid::IndexVector & >;

%warnfilter(401, 503) sgpp::combigrid::IndexVectorIterator;

// TODO(daissgr) Test in Python if the make_const_iterator can safely be
// removed without losing the iterator functionality. The problem with it: From
// what I can see in the swig source code this functionality was never inteded
// for python wrappers (mentions to %make_const_iterator only appear in the
// ruby context). Swig 4.1 and newer accordingly call it an unknown directive
// for the python wrappers (older swigs did not care but the statement probably
// never did work). Comment out for now to re-enable sgpp with newer swig
// versions and test them...
/* %make_const_iterator( IndexVectorIterator, IndexVectorRange ); */

%{
#include <sgpp/combigrid/tools/IndexVectorIterator.hpp>
#include <sgpp/combigrid/tools/IndexVectorRange.hpp>
%}

%include "combigrid/src/sgpp/combigrid/tools/IndexVectorIterator.hpp"
%include "combigrid/src/sgpp/combigrid/tools/IndexVectorRange.hpp"

%include "exception.i"
namespace sgpp {
  namespace combigrid {
    %exception IndexVectorIterator::__next__ {
      try
      {
        $action // calls %extend function next() below
      }
      catch (sgpp::combigrid::StopIterator)
      {
        PyErr_SetString(PyExc_StopIteration, "End of iterator");
        return NULL;
      }
    };

    %extend IndexVectorIterator
    {
      IndexVectorIterator*__iter__()
      {
        return $self;
      }

      IndexVector& __next__()
      {
        if (!($self->isAtEnd()))
        {
          // dereference the iterator and return reference to the object,
          // then increment
          return *($self->operator++());
        }
        throw sgpp::combigrid::StopIterator();
      }

      IndexVectorIterator* increment()
      {
        if (!($self->isAtEnd()))
        {
          // dereference the iterator and return reference to the object,
          // then increment
          $self->operator++();
          return $self;
        }
        throw sgpp::combigrid::StopIterator();
      }
    };

    %extend IndexVectorRange {
      IndexVectorIterator __iter__()
      {
        return $self->begin();
      }
    };
  }
}
