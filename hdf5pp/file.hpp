#ifndef HDF5_FILE_HPP
#define HDF5_FILE_HPP

#include "hdf5pp/handle.hpp"

#include <hdf5.h>
#include <stdexcept>

namespace hdf5 {
  class file : public handle {
  public:
    file() = default;
    
    explicit file(hid_t&& id) : handle(std::move(id)) { }
    
    //
    // static constructors
    //

    static file create(std::string const& name, unsigned int flags)
    {
      hid_t id = H5Fcreate(name.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
      if (id < 0) throw std::runtime_error("H5Fcreate failed");
      return file(std::move(id));
    }

    static file open(std::string const& name, unsigned int flags)
    {
      hid_t id = H5Fopen(name.c_str(), flags, H5P_DEFAULT);
      if (id < 0) throw std::runtime_error("H5Fopen failed");
      return file(std::move(id));
    }
   
    hsize_t get_filesize() const
    {
      hsize_t size;
      H5Fget_filesize(getId(), &size);
      return size;
    }
  };
}

#endif
