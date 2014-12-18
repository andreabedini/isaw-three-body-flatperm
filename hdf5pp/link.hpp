#ifndef LINK_HPP
#define LINK_HPP

#include "hdf5pp/handle.hpp"
#include <string>

namespace hdf5 {
  bool link_exists(handle const& loc,
		   std::string const& name,
		   hid_t lapl_id = H5P_DEFAULT)
  {
    htri_t r = H5Lexists(loc.getId(), name.c_str(), lapl_id);
    return (r > 0) ? true : false;
  }
}

#endif
