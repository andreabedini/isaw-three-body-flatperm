#ifndef DATATYPE_HPP
#define DATATYPE_HPP

#include <hdf5.h>

#include <cassert>
#include <string>

namespace hdf5 {
  //
  // forward declaration
  //
  template<typename T, typename Enable = void>
  struct datatype_from;

  //
  //
  //
  class datatype {
    hid_t id;

    //
    // friend classes
    //
    friend class dataset;

    template<typename T, typename Enable>
    friend struct datatype_from;

  public:
    operator hid_t() const { return id; }

    explicit datatype(hid_t&& id) : id(id)
    { }

    datatype(datatype const& o) : id(H5Tcopy(o.id))
    { };

    datatype(datatype&& o) : id(o.id)
    { o.id = 0; }

    ~datatype() { H5Tclose(id); }

    //////////////////////////////////////////////////////////////////////
    // static constructors
    //////////////////////////////////////////////////////////////////////

    static datatype create(H5T_class_t class_, size_t size)
    {
      hid_t id = H5Tcreate(class_, size);
      return datatype(std::move(id));
    }

    static datatype create_string()
    {
      hid_t id = H5Tcopy(H5T_C_S1);
      H5Tset_size(id, H5T_VARIABLE);
      return datatype(std::move(id));
    }

    static datatype create_string(size_t size)
    {
      hid_t id = H5Tcopy(H5T_C_S1);
      H5Tset_size(id, size);
      return datatype(std::move(id));
    }

    template<typename T>
    static datatype create_like(T const&)
    {
      return datatype_from<T>::value();
    }

    static datatype create_like(std::string const& t)
    {
      hid_t id = H5Tcopy(H5T_C_S1);
      H5Tset_size(id, t.size());
      return datatype(std::move(id));
    }

    //
    //
    //

    bool operator==(datatype const& rhs) const
    {
      htri_t res = H5Tequal(id, rhs.id);
      assert(res >= 0);
      return res;
    }

    //////////////////////////////////////////////////////////////////////
    //
    //////////////////////////////////////////////////////////////////////
  };
}

#endif
