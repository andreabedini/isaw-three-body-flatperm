#ifndef ATTRIBUTE_HPP
#define ATTRIBUTE_HPP

#include <hdf5pp/datatype.hpp>
#include <hdf5pp/dataspace.hpp>
#include <hdf5pp/handle.hpp>
#include <hdf5pp/utility/data.hpp>

#include <hdf5.h>

#include <stdexcept>

namespace hdf5 {
  class attribute : public handle {
  public:
    explicit attribute(hid_t&& id) : handle(std::move(id)) { }

    // static constructor

    static attribute create(handle const& loc, const char* name,
			    datatype const& type, dataspace const& space)
    {
      hid_t id = H5Acreate(loc.getId(), name, type, space,
			   H5P_DEFAULT, H5P_DEFAULT);
      assert(id >= 0);
      return attribute(std::move(id));
    }

    template<typename T>
    static attribute create(handle const& loc,
			    const char* name,
			    T const& t)
    {
      datatype type = datatype::create_like(t);
      dataspace space = dataspace::create_scalar();
      attribute attr = create(loc, name, type, space);
      attr.write(type, data(t));
      return attr;
    }

    void read(datatype const& mem_type, void * buf) {
      herr_t err = H5Aread(getId(), mem_type, buf);
      if (err < 0) throw std::runtime_error("H5Aread failed");
    }

    template<typename T>
    T read() {
      T t;
      read(datatype_from<T>::value(), data(t));
      return t;
    }

    void write(datatype const& mem_type, const void * buf) {
      herr_t err = H5Awrite(getId(), mem_type, buf);
      if (err < 0) throw std::runtime_error("H5Awrite failed");
    }

    template<typename T>
    void write(T const& t) {
      return write(datatype_from<T>::value(), data(t));
    }
  };

  //
  // Non member functions
  //
  
  attribute get_attribute(handle const& loc,
			  const char* obj_name,
			  const char* attr_name)
  {
    hid_t id = H5Aopen_by_name(loc.getId(), obj_name, attr_name,
			       H5P_DEFAULT, H5P_DEFAULT);
    if (id < 0) throw std::runtime_error("H5Aopen_by_name failed");
    return attribute(std::move(id));
  }

  attribute get_attribute(handle const& loc,
			  const char* attr_name)
  {
    return get_attribute(loc, ".", attr_name);
  }
}

#endif
