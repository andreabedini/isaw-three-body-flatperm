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
      // we need at least three points, i.e. two steps
      if (walk.size() < 2)
        return;

      const auto x = walk[walk.size() - 1];
      const auto origin = Walk::lattice_type::origin();

      // iterate over x's neighbours looking for new contacts
      for (auto const& y : Walk::lattice_type::get_neighbours(x)) {
        if (y != origin and walk.has_point(y) and not walk.has_bond(x, y)) {
          _m += 1;
        }
      }
    }

    template<typename Walk>
    void unregister_step(Walk const& walk)
    {
      if (walk.size() < 2)
        return;

      const auto x = walk[walk.size() - 1];
      const auto origin = Walk::lattice_type::origin();

      // iterate over x's neighbours
      for (auto const& y : Walk::lattice_type::get_neighbours(x)) {
        if (y != origin and walk.has_point(y) and not walk.has_bond(x, y))
          _m -= 1;
      }
    }

    int get() const { return _m; }
  };
}

#endif
