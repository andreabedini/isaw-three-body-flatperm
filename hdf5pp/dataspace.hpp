#ifndef HDF5_DATASPACE_HPP
#define HDF5_DATASPACE_HPP

#include <hdf5pp/utility/shape_of.hpp>
#include <hdf5pp/datatype.hpp>

#include <hdf5.h>

#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <array>
#include <cassert>
#include <initializer_list>

namespace hdf5 {
  using boost::enable_if;

  class dataspace {
    hid_t id;

  public:
    operator hid_t() const { return id; }

    explicit dataspace(hid_t&& id) : id(id)
    { }

    dataspace(dataspace const& o) : id(H5Scopy(o.id))
    { }

    dataspace(dataspace&& o) : id(o.id)
    { o.id = 0; }

    ~dataspace() { H5Sclose(id); }

    //////////////////////////////////////////////////////////////////////
    // static constructors
    //////////////////////////////////////////////////////////////////////    

    static dataspace
    create_scalar()
    {
      return dataspace(H5Screate(H5S_SCALAR));
    }

    template<typename T> static dataspace from(T const&);

    static dataspace
    create_simple(int rank, hsize_t const* current_dims, hsize_t const* maximum_dims = NULL)
    {
      return dataspace(H5Screate_simple(rank, current_dims, maximum_dims));
    }

    static dataspace
    create_simple(std::initializer_list<hsize_t> current_dims)
    {
      return create_simple(current_dims.size(), current_dims.begin());
    }

    template<typename T>
    static
    BOOST_CONCEPT_REQUIRES(((boost::Collection<T>)), (dataspace))
      create_simple(T const& current_dims)
    {
      return create_simple(current_dims.size(), current_dims.begin());
    }
    
    std::size_t get_simple_extent_ndims() const {
      return H5Sget_simple_extent_ndims(id);
    }

    hssize_t get_simple_extent_npoints() const {
      return H5Sget_simple_extent_npoints(id);
    }

    //
    // data selection
    //
    dataspace& select_all()
    {
      if (H5Sselect_all(id) < 0)
	throw 1;
      return *this;
    }

    dataspace& select_none()
    {
      if (H5Sselect_none(id) < 0 )
	throw 1;
      return *this;
    }

    bool selection_valid() const {
      // Returns a positive value, for TRUE, if the selection is
      // contained within the extent or 0 (zero), for FALSE, if it is
      // not. Returns a negative value on error conditions such as the
      // selection or extent not being defined.
      return H5Sselect_valid(id) > 0;
    }

    dataspace& select()
    { return *this; }
  };
}

#endif
