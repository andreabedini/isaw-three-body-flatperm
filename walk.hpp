/*
 * walk.hpp
 *
 */

#ifndef WALK_HPP__
#define WALK_HPP__

#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/next_prior.hpp>

#include <cassert>
#include <iostream>
#include <vector>

namespace mi = boost::multi_index;

namespace models {
template<typename Lattice>
class walk {
public:
  typedef Lattice lattice_type;
  typedef typename lattice_type::point point;

private:
  struct by_order;
  struct by_point;

  typedef boost::multi_index_container
  < point
  , mi::indexed_by
    < mi::random_access< mi::tag<by_order> >
    , mi::hashed_non_unique< mi::tag<by_point>, mi::identity<point>, boost::hash<point> >
    >
  > _structure_type;

  _structure_type _structure;

public:
  walk()
  {
    _structure.push_back(lattice_type::origin());
  }

  walk(unsigned int N)
  {
    _structure.reserve(N + 1);
    _structure.push_back(lattice_type::origin());
  }

  std::size_t size() const
  {
    return _structure.size() - 1;
  }

  point const& operator[](std::size_t i) const
  {
    return _structure[i];
  }

  decltype(_structure.begin()) begin() const { return _structure.begin(); }
  decltype(_structure.end())   end()   const { return _structure.end(); }

  decltype(_structure.rbegin()) rbegin() const { return _structure.rbegin(); }
  decltype(_structure.rend())   rend()   const { return _structure.rend();  }

  decltype(_structure.front()) front() const { return _structure.front(); }
  decltype(_structure.back())  back()  const { return _structure.back(); }

  //////////////////////////////////////////////////////////////////////
  
  int has_bond(point const& x, point const& y) const
  {
    auto range = _structure.template get<by_point>().equal_range(y);
    while (range.first != range.second) {
      auto i = _structure.template project<by_order>(range.first);
      if (i != _structure.begin() and x == *boost::prior(i))
    	return 1;
      if (x == *boost::next(i))
    	return 1;
      ++ range.first;
    }
    return 0;
  }

  int has_point(point const& x) const
  {
    return _structure.template get<by_point>().count(x);
  }

  void register_step(point const& x)
  {
    _structure.push_back(x);
  }

  void unregister_step()
  {
    assert(not _structure.empty());
    _structure.pop_back();
  }
};

template<typename Lattice>
std::ostream& operator<<(std::ostream& o, walk<Lattice> const& walk)
{
  for (auto i = walk.begin(); i != walk.end(); ) {
    o << *i;
    i++;
    if (i != walk.end())
      o << " -- ";

  }
  return o;
}
}

#endif // WALK_HPP__
