/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_TOOLS_SMARTPOINTER_HPP
#define SGPP_OPT_TOOLS_SMARTPOINTER_HPP

namespace sg {
  namespace opt {
    namespace tools {

      /**
       * Wrapper around a pointer which automatically frees memory on destruction.
       * It includes a counter, so that multiple SmartPointers can wrap the same pointer multiple times.
       * Only at the last SmartPointer destruction, the data is deleted.
       * (In C++11, std::unique_ptr could be used.)
       *
       * Adopted from http://www.codeguru.com/cpp/misc/samples/basicprogramming/article.php/c17495/
       * An-Efficient-Pointer-Wrapper-in-C-for-Scientific-Computation.htm.
       */
      template <class T>
      class SmartPointer {
        public:
          /// counter type
          typedef unsigned long ref_type;

          /**
           * Constructor.
           *
           * @param p     pointer to be wrapped
           */
          explicit SmartPointer(T* p = NULL) : pointer(p) {
            try {
              pref = new ref_type(1);
            } catch (...) {
              delete p;
              throw;
            }
          }

          /**
           * Copy constructor.
           * Increases the counter by one.
           */
          SmartPointer(const SmartPointer& r) throw() : pointer(r.pointer) {
            ++*r.pref;
            pref = r.pref;
          }

          /**
           * Assignment operator.
           * Calls dispose, copies the other SmartPointer and increases the counter by one.
           *
           * @param r     object to be assigned to
           */
          SmartPointer& operator=(const SmartPointer& r) throw() {
            if (this == &r) {
              return *this;
            }

            dispose();
            pointer = r.pointer;
            ++*r.pref;
            pref = r.pref;
            return *this;
          }

          /**
           * Assignment operator for pointers.
           * Calls dispose and essentially works like the constructor.
           */
          SmartPointer& operator=(T* p) {
            if (this->pointer == p) {
              return *this;
            }

            dispose();

            try {
              pref = new ref_type(1);
            } catch (...) {
              delete p;
              throw;
            }

            pointer = p;
            return *this;
          }

          /**
           * Destructor.
           * Calls dispose.
           */
          ~SmartPointer() {
            dispose();
          }

          /// SmartPointer with other type
          template <typename Y> friend class SmartPointer;

          /**
           * Constructor for other types.
           *
           * @param p     pointer to be wrapped
           */
          template <typename Y>
          SmartPointer(Y* py) {
            try {
              pref = new ref_type(1);
            } catch (...) {
              delete py;
              throw;
            }

            pointer = py;
          }

          /**
           * Copy constructor for other types.
           * Increases the counter by one.
           */
          template <typename Y>
          SmartPointer(const SmartPointer<Y>& r) {
            pointer = r.pointer;
            ++*r.pref;
            pref = r.pref;
          }

          /**
           * Assignment operator for other types.
           * Calls dispose, copies the other SmartPointer and increases the counter by one.
           *
           * @param r     object to be assigned to
           */
          template <typename Y>
          SmartPointer& operator=(const SmartPointer<Y>& r) {
            dispose();
            pointer = r.pointer;
            ++*r.pref;
            pref = r.pref;
            return *this;
          }

          /**
           * Assignment operator for pointers of other types.
           * Calls dispose and essentially works like the constructor.
           */
          template <typename Y>
          SmartPointer& operator=(Y* py) {
            if (this->pointer == py) {
              return *this;
            }

            dispose();

            try {
              pref = new ref_type(1);
            } catch (...) {
              delete py;
              throw;
            }

            pointer = py;
            return *this;
          }

          /**
           * Dispose old pointer and manage new pointer.
           *
           * @param p     new pointer to be wrapped
           */
          void reset(T* p = NULL) {
            if (pointer == p) {
              return;
            }

            if (--*pref == 0) {
              delete pointer;
            } else {
              try {
                pref = new ref_type;
              } catch (...) {
                ++*pref;
                delete p;
                throw;
              }
            }

            *pref = 1;
            pointer = p;
          }

          /**
           * @return  reference to dereferenced pointer
           */
          T& operator*() const throw() {
            return *pointer;
          }

          /**
           * @return  wrapped pointer
           */
          T* operator->() const throw() {
            return pointer;
          }

          /**
           * @return  wrapped pointer
           */
          T* get() const throw() {
            return pointer;
          }

          /**
           * @return  pointer counter
           */
          ref_type use_count() const throw() {
            return *pref;
          }

          /**
           * @return  true, if this is the only SmartPointer wrapping the pointer (counter equals 1)
           */
          bool unique() const throw() {
            return (*pref == 1);
          }

        private:
          /**
           * Decrements the counter and if it's zero, deletes pointer and counter.
           */
          void dispose() throw() {
            if (--*pref == 0) {
              delete pointer;
              delete pref;
            }
          }

          /// wrapped pointer
          T* pointer;
          /// pointer to counter
          ref_type* pref;
      };

    }
  }
}

#endif
