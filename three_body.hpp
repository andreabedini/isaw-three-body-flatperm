/*
 * three_body.hpp
 */

#ifndef THREE_BODY_HPP
#define THREE_BODY_HPP

#include <unordered_map>

namespace features {
  template<typename Lattice>
  class three_body {
    using point = typename Lattice::point;

    std::unordered_map<point, int, typename Lattice::hash> _faces;

    static std::pair<point, point>
    faces_from_step(point x, point y)
    {
      if (x[0] == y[0]) {
        // vertical step
        point m = x + y;
        point a = m + point{2,0};
        point b = m - point{2,0};
        return std::make_pair(a, b);
      } else {
        // horizontal step
        point a = (sum(x) % 2 == 0 ? x : y)*2 - point{0,1};
        point b = (sum(x) % 2 == 0 ? y : x)*2 + point{0,1};
        return std::make_pair(a, b);
      }
    }

  public:
    template<typename Walk>
    void register_step(Walk const& walk)
    {
      auto i = walk.rbegin();
      point y = *i++;              // last point
      point x = *i++;              // second last point

      // points a and b identify the faces adjacent to the step
      auto adjacent_faces = faces_from_step(x, y);
      _faces[adjacent_faces.first] ++;
      _faces[adjacent_faces.second] ++;
    }

    template<typename Walk>
    void unregister_step(Walk const& walk)
    {
      auto i = walk.rbegin();
      point y = *i++;             // last point
      point x = *i++;             // second last point

      auto adjacent_faces = faces_from_step(x, y);
      _faces[adjacent_faces.first] --;
      _faces[adjacent_faces.second] --;
    }

    int get() const {
      // count the number of faces we have seen, each with multiplity one
      int c = 0;
      for (auto const& el : _faces) {
        if (el.second > 0)
          c ++;
      }
      return c;
    }
  };
}

#endif // THREE_BODY_HPP

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=2 ts=2 : */
