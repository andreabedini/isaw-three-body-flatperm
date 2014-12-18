#ifndef HDF5_DATASET_HPP
#define HDF5_DATASET_HPP

#include "hdf5pp/handle.hpp"
#include <hdf5pp/dataspace.hpp>
#include <hdf5pp/datatype.hpp>
#include <hdf5pp/link.hpp>
#include <hdf5pp/utility/data.hpp>

#include <hdf5.h>

#include <stdexcept>

namespace hdf5 {
  template<typename Base>
  struct __link_create {
    hid_t lcpl_id;
  public:
    __link_create() : lcpl_id(H5Pcreate(H5P_LINK_CREATE)) {}
    ~__link_create() { H5Pclose(lcpl_id); }
    hid_t get_link_create() const { return lcpl_id; };
  };

  template<typename Base>
  struct __dataset_create {
    hid_t dcpl_id;
    __dataset_create() : dcpl_id(H5Pcreate(H5P_DATASET_CREATE)) {}
    ~__dataset_create() { H5Pclose(dcpl_id); }

    template<typename Collection>
    Base& set_chunk(Collection const& dims) {
      H5Pset_chunk(dcpl_id, dims.size(), dims.begin());
      return *static_cast<Base*>(this);
    }
  };

  template<typename Base>
  struct __dataset_access {
    hid_t dapl_id;
    __dataset_access() : dapl_id(H5Pcreate(H5P_DATASET_ACCESS)) {}
    ~__dataset_access() { H5Pclose(dapl_id); }
  };

  // forward declation
  class dataspace;

  class dataset : public handle {
  public:
    dataset() = default;    

    explicit dataset(hid_t&& id) : handle(std::move(id)) { }

    // dataset(dataset&& o) : handle(std::move(o)) { }

    //////////////////////////////////////////////////////////////////////
    // H5Dcreate
    //////////////////////////////////////////////////////////////////////

    class create
      : public __link_create<create>
      , public __dataset_create<create>
      , public __dataset_access<create>
    {
      const handle& loc;
      const char* name;
      const datatype type;
      const dataspace space; 
      friend class dataset;
   public:
      create(handle const& loc, const char* name,
	     datatype const& type, dataspace const& space)
	: loc(loc), name(name), type(type), space(space)
      {}
    };
    
    dataset(create const& params)
      : handle(H5Dcreate(params.loc.getId(), params.name,
			 params.type, params.space,
			 params.lcpl_id, params.dcpl_id, params.dapl_id))
    { }

    //////////////////////////////////////////////////////////////////////
    // H5Dopen
    //////////////////////////////////////////////////////////////////////

    static dataset open(handle const& loc, const char* name)
    {
      return dataset(H5Dopen(loc.getId(), name, H5P_DEFAULT));
    }

    //////////////////////////////////////////////////////////////////////
    //
    //////////////////////////////////////////////////////////////////////

    datatype get_type() const
    {
      return datatype(H5Dget_type(getId()));
    }

    dataspace get_space() const
    {
      return dataspace(H5Dget_space(getId()));
    }

    hsize_t get_storage_size() const
    {
      return H5Dget_storage_size(getId());
    }

    //////////////////////////////////////////////////////////////////////
    // transfer operations
    //////////////////////////////////////////////////////////////////////

    herr_t read(datatype const& mem_type, dataspace const& mem_space, void *buf)
    {
      return H5Dread(getId(), mem_type, mem_space, H5S_ALL, H5P_DEFAULT, buf);
    }

    template<typename T>
    const T read()
    {
      datatype type   = datatype_from<T>::value();
      dataspace space = dataspace::create_scalar();

      T t;

      read(type, space, data(t));

      return t;
    }

    void write(datatype const& mem_type, dataspace const& mem_space,
	       const void *buf)
    {
      H5Dwrite(getId(), mem_type, mem_space, H5S_ALL, H5P_DEFAULT, buf);
    }

    template<typename T>
    void write(dataspace const& mem_space, T const* buf)
    {
      H5Dwrite(getId(), datatype_from<T>::value(), mem_space,
	       H5S_ALL, H5P_DEFAULT, buf);
    }

    const std::string read_string()
    {
      datatype type   = datatype::create_string();
      dataspace space = dataspace::create_scalar();

      char *p[1];

      read(type, space, p);

      std::string s(p[0]);

      if (p[0]) free(p[0]);

      return s;
    }

    void write_string(std::string const& s)
    {
      datatype type   = hdf5::datatype::create_string();
      dataspace space = hdf5::dataspace::create_scalar();
      void const* p[1] = { data(s) };
      write(type, space, p);
    }
  };
}

#endif
