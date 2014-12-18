/*
 * lattice.hpp
 *
 */

#ifndef LATTICE_HPP__
#define LATTICE_HPP__

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

namespace lattices {
  template<unsigned int N, typename T = int>
  struct point : public std::array<T, N> {
    template<typename U>
    struct rebind { typedef point<N, U> result_type; };

    using std::array<T, N>::begin;
    using std::array<T, N>::end;

    point() {
      std::fill(begin(), end(), 0);
    }

    template<typename U>
    point(std::initializer_list<U> list) {
      std::copy_n(list.begin(), N, begin());
    }

    template<typename U>
    point(const point<N, U>& other) {
      std::copy_n(other.begin(), N, begin());
    }
    
    point& operator+=(const point& rhs) {
      std::transform(begin(), end(), rhs.begin(), begin(),
        std::plus<T>());
      return *this;
    }

    const point operator+(const point& rhs) const {
      point result(*this);
      result += rhs;
      return result;
    }

    point& operator-=(const point& rhs) {
      std::transform(begin(), end(), rhs.begin(), begin(),
        std::minus<T>());
      return *this;
    }

    const point operator-(const point& rhs) const {
      point result(*this);
      result -= rhs;
      return result;
    }

    point& operator*=(int k) {
      std::transform(begin(), end(), begin(),
          [k](T x) { return x * k; });
      return *this;
    }

    const point operator*(int k) const {
      point result(*this);
      result *= k;
      return result;
    }
 
    point& operator/=(int k) {
      std::transform(begin(), end(), begin(),
          [k](T x) { return x / k; });
      return *this;
    }

    const point operator/(int k) const {
      point result(*this);
      result /= k;
      return result;
    }
  };
  
  template<unsigned int N, typename T>
  std::ostream& operator<<(std::ostream& o, const point<N, T>& rhs) {
    auto i = std::begin(rhs);
    o << "(";
    while (true) {
      o << *i++;
      if (i == std::end(rhs))
        break;
      else
        o << ",";
    }
    return o << ")";
  }

  template<unsigned int N, typename T>
  T sum(const point<N, T>& p)
  {
    return std::accumulate(p.begin(), p.end(), T(0));
  }

  template<unsigned int N, typename T>
  T norm_square(const point<N, T>& p)
  {
    return std::inner_product(p.begin(), p.end(), p.begin(), T(0));
  }

  template<typename T>
  T cross(const point<2, T>& lhs, const point<2, T>& rhs)
  {
    return lhs[0] * rhs[1] - lhs[1] * rhs[0];
  }

  //////////////////////////////////////////////////////////////////////

  struct hexagonal {
    static const unsigned int dimensionality = 2;
    static const unsigned int coordination = 3;

    struct point : public lattices::point<2> {
      using lattices::point<2>::point;
    };

    lattices::point<2> to_xy(point const& p) {
      // FIXME
      return p; 
    }

    int norm_square(point const& p) {
      return norm_square(to_xy(p));
    }

    struct hash {
      std::size_t operator()(point const& p) const {
        std::size_t const h1 = std::hash<int>()(p[0]);
        std::size_t const h2 = std::hash<int>()(p[1]);
        return h1 ^ (h2 << 1);
      }
    };

    static std::vector<point> get_neighbours(point const& c)
    {
      if (sum(c) % 2 == 0)
        return { c + point{1,0}, c + point{-1,0}, c + point{0,1} };
      else
        return { c + point{1,0}, c + point{-1,0}, c + point{0,-1} };
    }

    static const point origin() { return point{0, 0}; };
  };
}

#endif // LATTICE_HPP__

// vim: noai:ts=2:sw=2

