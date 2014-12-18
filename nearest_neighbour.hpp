/*
 * $Id: nearest_neighbour.hpp,v 1.1 2012/01/11 07:32:53 abedini Exp $
 *
 * $Log: nearest_neighbour.hpp,v $
 * Revision 1.1  2012/01/11 07:32:53  abedini
 * Initial revision
 *
 *
 */

#ifndef NEAREST_NEIGHBOR_HPP
#define NEAREST_NEIGHBOR_HPP

namespace features {
  class nearest_neighbour {
    int _m = 0;

  public:
    template<typename Walk>
    void register_step(Walk const& walk)
    {
      // fetch the last point
      auto const& x = walk.back();

      // iterate over x's neighbours looking for new contacts
      for (auto const& y : Walk::lattice_type::get_neighbours(x)) {
        if (walk.has_point(y) and not walk.has_bond(x, y))
          _m += 1;
      }
    }

    template<typename Walk>
    void unregister_step(Walk const& walk)
    {
      // fetch the last point
      auto const& x = walk.back();

      // iterate over x's neighbours
      for (auto const& y : Walk::lattice_type::get_neighbours(x)) {
        if (walk.has_point(y) and not walk.has_bond(x, y))
          _m -= 1;
      }
    }

    int get() const { return _m; }
  };
}

#endif
