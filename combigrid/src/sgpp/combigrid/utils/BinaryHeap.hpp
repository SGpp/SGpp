// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Implements a max heap using a binary heap. We need this as a priority queue implementation with a
 * updatePriority() function because we do not want to have boost as a dependency.
 */
template <typename T, typename Comparator = std::less<T>>
class BinaryHeap {
 public:
  struct HeapElement {
    T elem;
    size_t idx;

    HeapElement(T e, size_t i) : elem(e), idx(i) {}
  };

  typedef std::shared_ptr<HeapElement> HeapElementPointer;

  struct Handle {
    HeapElementPointer elementPointer;

    explicit Handle(HeapElementPointer elementPointer) : elementPointer(elementPointer) {}
    T& operator*() { return elementPointer->elem; }
  };

  typedef Handle handle_type;

 private:
  std::vector<HeapElementPointer> data;

  size_t get_lchild_index(size_t nodeIndex) { return 2 * nodeIndex + 1; }

  size_t get_rchild_index(size_t nodeIndex) { return 2 * nodeIndex + 2; }

  size_t get_parent_index(size_t nodeIndex) { return (nodeIndex - 1) / 2; }

  bool isSmaller(T const& first, T const& second) { return Comparator()(first, second); }
  bool isGreater(T const& first, T const& second) { return isSmaller(second, first); }

  inline void swapElements(HeapElementPointer& a, HeapElementPointer& b) {
    /*T tmp = a->elem;
    a->elem = b->elem;
    b->elem = tmp;*/
    std::swap(data[a->idx], data[b->idx]);
    std::swap(a->idx, b->idx);
  }

  HeapElementPointer shift_up(size_t node_index) {
    if (node_index != 0) {
      size_t parent_index = get_parent_index(node_index);
      if (isGreater(data[node_index]->elem, data[parent_index]->elem)) {
        swapElements(data[node_index], data[parent_index]);
        return shift_up(parent_index);
      }
    }
    return data[node_index];
  }

  HeapElementPointer shift_down(size_t node_index) {
    size_t lchild_index = get_lchild_index(node_index);
    size_t rchild_index = get_rchild_index(node_index);
    size_t max_index;
    // find the max value from two children
    if (rchild_index >= data.size()) {
      if (lchild_index >= data.size())
        return data[node_index];
      else
        max_index = lchild_index;
    } else {
      if (!isSmaller(data[lchild_index]->elem, data[rchild_index]->elem))
        max_index = lchild_index;
      else
        max_index = rchild_index;
    }
    // compare node value with max value
    if (isSmaller(data[node_index]->elem, data[max_index]->elem)) {
      swapElements(data[node_index], data[max_index]);
      return shift_down(max_index);
    }
    return data[node_index];
  }

 public:
  ~BinaryHeap() { data.clear(); }

  BinaryHeap() {}

  explicit BinaryHeap(size_t size) { data.reserve(size); }

  explicit BinaryHeap(std::vector<T> v) {
    data.reserve(v.size());
    for (auto it = v.begin(); it != v.end(); ++it) {
      push(*it);
    }
  }

  bool empty() { return (data.size() == 0); }

  HeapElementPointer get(size_t index) { return data[index]; }

  void set(size_t index, T value) { data[index]->elem = value; }

  void clear() { data.clear(); }

  T const& top() {
    if (empty()) throw std::string("Heap is empty");
    return data[0]->elem;
  }

  Handle push(T elem) {
    data.push_back(HeapElementPointer(new HeapElement(elem, data.size())));
    return Handle(shift_up(data.size() - 1));
  }

  void pop() {
    if (empty()) {
      throw std::string("Heap is empty");
    }
    swapElements(data[0], data[data.size() - 1]);
    data.pop_back();
    if (data.size() > 1) {
      shift_down(0);
    }
  }

  void update(Handle h_in, T const& value) {
    *h_in = value;
    HeapElementPointer h_tmp = shift_up(h_in.elementPointer->idx);
    shift_down(h_tmp->idx);
  }

  void update(Handle h_in) {
    HeapElementPointer h_tmp = shift_up(h_in->idx);
    shift_down(h_tmp->idx);
  }

  // for debug only
  void print() {
    for (auto it = data.begin(); it != data.end(); ++it)
      std::cout << "[" << (*it)->elem << ", " << (*it)->idx << "]" << std::endl;
    std::cout << std::endl;
  }
};
} /* namespace combigrid */
} /* namespace sgpp */
