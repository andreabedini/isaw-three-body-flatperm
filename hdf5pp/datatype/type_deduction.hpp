#ifndef TYPE_DEDUCTION_HPP
#define TYPE_DEDUCTION_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/type_traits/remove_all_extents.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include "hdf5pp/datatype.hpp"
#include "hdf5pp/utility/shape_of.hpp"

namespace hdf5 {
  template<typename T>
  struct datatype_from<T, typename boost::enable_if<boost::is_array<T>>::type>
  {
    typedef typename boost::remove_all_extents<T>::type base_type;
    static datatype value()
    {
      auto const shape = utility::Shape<T>::as_array();
      return H5Tarray_create(datatype_from<base_type>::value(),
			     shape.size(), shape.data());
    }
  };

  template<typename T>
  struct datatype_from<T, typename boost::enable_if<boost::is_reference<T> >::type>
  {
    static datatype value()
    {
      return datatype_from<typename boost::remove_reference<T>::type>::value();
    }
  };

  template<typename T>
  struct datatype_from<T, typename boost::enable_if<boost::is_const<T> >::type>
  {
    static datatype value()
    {
      return datatype_from<typename boost::remove_const<T>::type>::value();
    }
  };
}

#endif
