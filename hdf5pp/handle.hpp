#ifndef HANDLE_HPP
#define HANDLE_HPP

#include <hdf5.h>

namespace hdf5 {
  class handle {
    hid_t _id;

  public:
    handle() : _id(0) { }

    explicit handle(hid_t&& id) : _id(id) { }
    
    handle(handle const& o) : _id(o._id) {
      if (_id) H5Iinc_ref(_id);
    }
    
    handle(handle&& o) : _id(std::move(o._id)) { }

    virtual ~handle() {
      if (_id) H5Idec_ref(_id);
    }
    
    handle& operator=(handle const& o)
    {
      if (_id) H5Idec_ref(_id);
      _id = o._id;
      if (_id) H5Iinc_ref(_id);
      return *this;
    }

    handle& operator=(handle&& o)
    {
      using namespace std;
      swap(_id, o._id);
      return *this;
    }

    bool is_valid() const { return H5Iis_valid(_id) == true; }
    hid_t getId() const { return _id; }

    void flush()
    {
      herr_t err = H5Fflush(getId(), H5F_SCOPE_GLOBAL);
      if (err < 0) throw std::runtime_error("H5Fflush failed");      
    }
  };
}

#endif
