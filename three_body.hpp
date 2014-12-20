/*
 * three_body.hpp
 */

#ifndef THREE_BODY_HPP
#define THREE_BODY_HPP

#include <unordered_map>
#include <boost/container/static_vector.hpp>

namespace features {
  template<typename Lattice>
  class three_body {
    std::array<int, Lattice::coordination> k;

    struct __data {
      int contacts;
      boost::container::static_vector<int, 6> visits;
    };

    using point = typename Lattice::point;
    std::unordered_map<point, __data, typename Lattice::hash> _faces;

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

    void do_face_in(__data& data, int n) {
      if (data.contacts == 0 or (data.visits.back() < n - 1)) {
        data.contacts ++;
        -- k[data.contacts-1];
        ++ k[data.contacts];
      }

      data.visits.push_back(n);
    }

    void do_face_out(__data& data, int n) {
      data.visits.pop_back();

      // check if we are the last visit
      if (data.visits.empty() or data.visits.back() < n - 1) {
        -- k[data.contacts];
        ++ k[data.contacts-1];
        data.contacts --;
      }

      assert(data.contacts >= 0);

      if (n == 1) {
        assert(data.contacts == 0);
        assert(data.visits.empty());
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
      const auto adjacent_faces = faces_from_step(x, y);
      const auto n = walk.size();

      do_face_in(_faces[adjacent_faces.first], n);
      do_face_in(_faces[adjacent_faces.second], n);
    }

    template<typename Walk>
    void unregister_step(Walk const& walk)
    {
      auto i = walk.rbegin();
      point y = *i++;             // last point
      point x = *i++;             // second last point

      const auto adjacent_faces = faces_from_step(x, y);
      const auto n = walk.size();

      do_face_out(_faces[adjacent_faces.first], n);
      do_face_out(_faces[adjacent_faces.second], n);

      if (walk.size() == 1) {
        assert(get_N2() == 0);
        assert(get_N3() == 0);
      }
    }

    int get_N2() const {
      return k[2];
    }

    int get_N3() const {
      return k[3];
    }
  };
}

#endif // THREE_BODY_HPP

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=2 ts=2 : */
