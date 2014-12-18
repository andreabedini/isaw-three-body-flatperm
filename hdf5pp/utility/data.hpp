#ifndef DATA_HPP
#define DATA_HPP

#include <string>

namespace hdf5 {
  template<typename T>
  void* data(T& t) { return &t; }

  template<typename T>
  void const* data(T const& t) { return &t; }

  template<>
  void const* data(std::string const& s) { return s.data(); }
}

#endif
