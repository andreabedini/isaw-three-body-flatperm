#ifndef SHAPE_OF_HPP
#define SHAPE_OF_HPP

#include <array>

namespace hdf5 {
  namespace utility {
    //////////////////////////////////////////////////////////////////////
    
    template<int... Is> struct SizeList { 
      typedef std::array<size_t, sizeof...(Is)> array_type;
      static array_type as_array() { return array_type{Is...} ; }
    };
    
    template<typename SizeList, int I>
    struct Insert;
    
    template<int... Is, int I>
    struct Insert<SizeList<Is...>, I> {
      typedef SizeList<I, Is...> result;
    };
    
    //////////////////////////////////////////////////
    
    template<typename T>
    struct Shape {
      typedef SizeList<> result;
    };
    
    template<typename T, int N>
    struct Shape<T[N]> {
      typedef typename Insert<typename Shape<T>::result, N>::result result;
    };
    
    //////////////////////////////////////////////////
    
    template<typename T>
    auto shape_of(T const&) -> typename Shape<T>::result::array_type
    {
      return Shape<T>::result::as_array();
    }
    
    //////////////////////////////////////////////////////////////////////
  }
}

#endif
