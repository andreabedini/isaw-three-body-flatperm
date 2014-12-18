#ifndef MY_ARRAY_HPP
#define MY_ARRAY_HPP

#include "hdf5pp/hdf5.hpp"

#include <algorithm>
#include <vector>

template<typename ValueType, size_t D>
class my_array {
public:
  static const size_t NumDims = D;

  typedef size_t size_type;

  typedef ValueType        value_type;
  typedef ValueType&       reference;
  typedef ValueType const& const_reference;

  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

private:
  size_type __extents[NumDims];
  std::vector<value_type> __data;

  template<size_t I>
  struct helper {
    static_assert(NumDims >= I, "more indices than dimensions");

    my_array& array;
    const size_t offset;

    operator reference() { return array.__data[offset]; }

    template<typename T>
    reference operator=(T const& t) {
      return array.__data[offset] = t;
    }

    helper<I + 1> operator[](size_t i) const {
      return helper<I + 1>{array, i + array.__extents[I] * offset};
    }
  };

public:
  my_array()
  {
    std::fill(std::begin(__extents), std::end(__extents), 0);
  }

  my_array(std::initializer_list<unsigned int> extents)
  {
    std::copy_n(std::begin(extents), NumDims, std::begin(__extents));
    __data.resize(num_elements());
  }

  template<typename ExtentList>
  my_array(ExtentList const& extents)
  {
    std::copy_n(std::begin(extents), NumDims, std::begin(__extents));
    __data.resize(num_elements());
  }

  size_t num_dimensions() const
  {
    return NumDims;
  }

  size_t num_elements() const
  {
    return std::accumulate(std::begin(__extents), std::end(__extents),
			   1, std::multiplies<size_type>());
  }

  iterator begin() { return __data.begin(); }
  iterator end()   { return __data.end(); }

  const_iterator begin() const { return __data.begin(); }
  const_iterator end()   const { return __data.end(); }

  const value_type* data() const { return &__data[0]; }
  value_type*       data()       { return &__data[0]; }

  const size_type*
  shape() const { return &__extents[0]; }

  helper<1> operator[](size_type i)
  {
    return helper<0>{*this, 0}[i];
  }

  template<typename IndexList>
  value_type& operator()(IndexList const& indices)
  {
    size_type offset = indices[0];
    for (size_type n = 1; n != NumDims; ++n) {
      assert(indices[n] < __extents[n]);
      offset = indices[n] + __extents[n] * offset;
    }

    return __data[offset];
  }
};

//////////////////////////////////////////////////////////////////////

namespace hdf5 {
  template <typename ValueType, size_t NumDims>
  void load(handle loc, my_array<ValueType, NumDims>& h,
	    const char* name)
  {
    datatype type = datatype_from<ValueType>::value();
    std::array<hsize_t, NumDims> extents;
    std::copy(h.shape(), h.shape() + NumDims, extents.begin());
    dataspace space = dataspace::create_simple(extents);

    ( link_exists(loc, name)
      ? dataset::open(loc, name)
      : dataset::create(loc, name, type, space) )
      .read(type, space, h.data());
  }

  template <typename ValueType, size_t NumDims>
  void save(handle loc, my_array<ValueType, NumDims> const& h,
	    const char* name)
  {
    datatype type = datatype_from<ValueType>::value();
    std::array<hsize_t, NumDims> extents;
    std::copy(h.shape(), h.shape() + NumDims, extents.begin());
    dataspace space = dataspace::create_simple(extents);

    ( link_exists(loc, name)
      ? dataset::open(loc, name)
      : dataset::create(loc, name, type, space) )
      .write(type, space, h.data());
  }
}

#endif
