#ifndef HDF5_LITE_HPP
#define HDF5_LITE_HPP

#include <hdf5_hl.h>
#include <cassert>

namespace hdf5 {
  const char *datatype_to_text(datatype const& type)
  { 
    char *buf = NULL;
    size_t size = 0;
    hid_t id;
    id = H5LTdtype_to_text(type._id, NULL, H5LT_DDL, &size);
    assert(id >= 0);
    buf = new char[size];
    id = H5LTdtype_to_text(type._id, buf, H5LT_DDL, &size);
    assert(id >= 0);
    return buf;
 }
}

#endif
