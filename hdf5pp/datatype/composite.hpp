#ifndef COMPOSITE_HPP
#define COMPOSITE_HPP

#include <hdf5pp/datatype/datatype.hpp>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>

#define HDF5_NAME_(type) BOOST_PP_CAT(__, type)

#define HDF5_COMPOSITE_INSERT(r, type, field)				\
  H5Tinsert(HDF5_NAME_(type), BOOST_PP_STRINGIZE(field),		\
	    offsetof(type, field),					\
	    hdf5::datatype_from<decltype(type().field)>::value());

#define HDF5_COMPOSITE_DEFINE(type, seq)				\
  hid_t HDF5_NAME_(type) = H5Tcreate(H5T_COMPOUND, sizeof(type));	\
  BOOST_PP_SEQ_FOR_EACH(HDF5_COMPOSITE_INSERT, type, seq);		\

#define HDF5_COMPOSITE(type)			\
  hdf5::datatype(std::move(HDF5_NAME_(type)))

#endif
