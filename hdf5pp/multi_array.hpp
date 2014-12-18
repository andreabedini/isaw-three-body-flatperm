#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

#include "hdf5pp/handle.hpp"

#include <array>
#include <boost/multi_array.hpp>

namespace hdf5 {
  template <typename ValueType, std::size_t NumDims, typename Allocator>
  void load(handle const& loc,
	    boost::multi_array<ValueType, NumDims, Allocator>& h,
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
  
  template <typename ValueType, std::size_t NumDims, typename Allocator>
  void save(handle const& loc,
	    boost::multi_array<ValueType, NumDims, Allocator> const& h,
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
