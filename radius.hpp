/*
 * radius.hpp
 * 
 */

#ifndef RADIUS_HPP
#define RADIUS_HPP

namespace features {
  template<typename point>
  struct radius {
    typename point::template rebind<int64_t>::result_type B;
    int64_t C = 0;
    
    template<typename Walk>
    void register_step(Walk const& walk) {
      point p = walk.back();
      
      B += p;
      C += norm_square(p);
    }

    template<typename Walk>
    void unregister_step(Walk const& walk) {
      point p = walk.back();

      B -= p;
      C -= norm_square(p);
    }

    int64_t get_CM_norm_square() const {
      return norm_square(B);
    }

    int64_t get_norm_square_sum() const {
      return C;
    }

  };
}

#endif
