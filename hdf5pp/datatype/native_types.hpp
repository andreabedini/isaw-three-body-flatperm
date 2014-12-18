#ifndef NATIVE_TYPES_HPP
#define NATIVE_TYPES_HPP

#include <string>

#include "hdf5pp/datatype.hpp"
#include "hdf5pp/datatype/type_deduction.hpp"

namespace hdf5 {

#define __HDF5_MAP_TYPE( type, hdf5_type )\
  template<> struct datatype_from<type>\
  { static datatype value() { return datatype(H5Tcopy(hdf5_type)); }}

  __HDF5_MAP_TYPE(int, H5T_NATIVE_INT);
  __HDF5_MAP_TYPE(unsigned int, H5T_NATIVE_UINT);

  __HDF5_MAP_TYPE(long int, H5T_NATIVE_LONG);
  __HDF5_MAP_TYPE(unsigned long int, H5T_NATIVE_ULONG);

  __HDF5_MAP_TYPE(long long, H5T_NATIVE_LLONG);
  __HDF5_MAP_TYPE(unsigned long long, H5T_NATIVE_ULLONG);

  __HDF5_MAP_TYPE(float, H5T_NATIVE_FLOAT);
  __HDF5_MAP_TYPE(double, H5T_NATIVE_DOUBLE);
  __HDF5_MAP_TYPE(long double, H5T_NATIVE_LDOUBLE);

#undef __HDF5_MAP_TYPE
}

#endif
