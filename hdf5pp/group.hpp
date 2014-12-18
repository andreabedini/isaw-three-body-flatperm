#ifndef HDF5_GROUP_HPP
#define HDF5_GROUP_HPP

#include "hdf5pp/handle.hpp"

#include <hdf5.h>

namespace hdf5 {
  class group : public handle {
  public:
    group() = default;

    explicit group(hid_t&& id) : handle(std::move(id)) { }

    //////////////////////////////////////////////////////////////////////
    // static constructors
    //////////////////////////////////////////////////////////////////////    
    
    static
    group open(hid_t loc, std::string const& name)
    {
      return H5Gopen(loc, name.c_str(), H5P_DEFAULT);
    }
  };
}

#endif
