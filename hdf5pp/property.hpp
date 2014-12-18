#ifndef HDF5_PROPERTY_HPP
#define HDF5_PROPERTY_HPP

#include "hdf5pp/datatype.hpp"

#include <hdf5.h>

#include <initializer_list>

namespace hdf5 {
  struct dataset_creation_property_list {
    hid_t dcpl_id;

    dataset_creation_property_list()
      : dcpl_id(H5Pcreate(H5P_DATASET_CREATE))
    { }

    template<typename X>
    dataset_creation_property_list&
    set_fill_value(X const& value)
    {
      H5Pset_fill_value(dcpl_id, datatype::from(value), (const void *) &value);
      return *this;
    }
    
    template<typename... Args>
    dataset_creation_property_list&
    set_chunk(Args... args)
    {
      hsize_t dims[] = { args... };
      H5Pset_chunk(dcpl_id, sizeof...(Args), dims);
      return *this;
    }
  };
}

#endif
