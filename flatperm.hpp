/*
 * flatperm.hpp
 *
 */

#ifndef FLATPERM_HPP
#define FLATPERM_HPP

#include "hdf5pp/hdf5.hpp"

#include <array>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <initializer_list>
#include <random>

#include "my_array.hpp"

#include <boost/thread/thread.hpp>

template<typename T, std::size_t N>
struct better_array : public std::array<T, N> {
  using base_type = std::array<T, N>;
  // initialize with zeroes
  better_array() {
    base_type::fill(T(0));
  }
  // initialize with initialization_list syntax
  better_array(std::initializer_list<T> a) {
    std::copy_n(std::begin(a), N, base_type::begin());
  }
};

namespace algorithm {
template<unsigned int D, typename RandomGenerator>
struct flatperm {
  //////////////////////////////////////////////////
  // random number generator and distributions
  //////////////////////////////////////////////////
  RandomGenerator& rng;

  //////////////////////////////////////////////////
  // flatperm histograms
  //////////////////////////////////////////////////

  using indices_type = better_array<unsigned int, D>;

  indices_type indices;
  const indices_type extents;

  my_array<long double, D> sW, Se;
  my_array<uint64_t, D> Sn, Enr, Pru;

  const long double mu;

  flatperm(std::initializer_list<unsigned int> extents_, double mu, RandomGenerator& rng)
    : rng(rng)
    , extents(extents_)
    , sW(extents), Se(extents), Sn(extents)
    , Enr(extents), Pru(extents)
    , mu(mu)
  {
    std::cerr << "Flatperm initialized, ";
    std::cerr << "extents ";
    for (auto x : extents)
      std::cerr << x << " ";
    std::cerr << ", weight renormalization set to " << mu << ".\n";
  }

  ////////////////////////////////////////////////////
  //
  ////////////////////////////////////////////////////
  template<typename T>
  void run(T* instance, unsigned int Snew)
  {
    std::uniform_real_distribution<double> uniform01;

    auto const Nmax = extents[0] - 1;
    auto const& walk_size = indices[0];

    //////////////////////////////////////////////////
    // history management
    //////////////////////////////////////////////////

    typedef decltype(instance->atmosphere()) atmosphere_type;

    struct mark {
      long unsigned int n;
      long double W;
      atmosphere_type enrichments;
    };

    std::vector<mark> history;
    history.reserve(extents[0]);

    auto last_enrichment = [&] {
      return history.empty() ? 0 : history.back().n;
    };

    const double delay = 0.1;

    uint64_t S = sW.data()[0];
    uint64_t Smax = S + Snew;

    std::cerr << "I already have " << S << " tours, starting " << Snew
	      << " up to " << Smax << "\n";

    while (S < Smax) {
      assert(history.empty());
      assert(walk_size == 0);

      // begin new tour
      indices.fill(0);

      S += 1;
      long double W = 1;

      sW(indices) += W;
      Sn(indices) += 1;
      Se(indices) += 1;

      while (true) {
	boost::this_thread::interruption_point();

	// Step 1 - get the atmosphere
	auto atmo = instance->atmosphere();

	// Step 2 - prune or enrich
	// The following piece compute 'copies' and possibily updates 'W'

	size_t copies = 0;

	if (walk_size < Nmax and not atmo.empty() and delay * walk_size < S) {
	  long double const Srel = S - std::floor(delay * walk_size);
	  long double const target_weight = sW(indices) / Srel;
	  long double const tw_correction = Se(indices) / Srel;
	  long double const ratio = W / target_weight / tw_correction;

	  if (ratio < 1.0) {
	    // probabilistic pruning
	    if (uniform01(rng) < ratio) {
	      copies = 1;
	      W /= ratio;
	    } else {
	      copies = 0;
	      W = 0;
	    }
	  } else {
	    copies = std::min(atmo.size(), (size_t) floor(ratio));
	    W /= copies;
	  }
	} else {
	  copies = 0;
	  W = 0;
	}

  	// Standard RR step. It is alright doing it here, after having
	// determined how many copies but before the actual
	// pruning/enrichment, basically because the weight has not to
	// change between "history" operations (push/pull).
	W *= atmo.size() / mu;

	// Step 3 - shrink and reload (if needed)
	if (copies == 0) {
	  // stats
	  Pru(indices) ++;

	  // Shrink the walk
	  auto const n_last = last_enrichment();
	  while (walk_size > n_last)
	    instance->unregister_step();

	  // check if we finished a tour
	  if (history.empty())
	    break;

	  W = history.back().W;
	} else {
	  assert(copies > 0);
	  assert(not atmo.empty());

	  // stats
	  Enr(indices) += copies - 1;

	  // sample 'copies' from the atmosphere
	  shuffle(begin(atmo), end(atmo), rng);
	  atmosphere_type enrichments;
	  enrichments.reserve(copies);
	  copy_n(begin(atmo), copies, back_inserter(enrichments));

	  history.push_back(mark{walk_size, W, std::move(enrichments)});
	}

	assert(not history.empty());
	assert(not history.back().enrichments.empty());

	auto const next_point = history.back().enrichments.back();
	history.back().enrichments.pop_back();

	// Step 4c - clean the history
	if (history.back().enrichments.empty())
	  history.pop_back();

	// Step 4 - Add a new node
	instance->register_step(next_point, W);

	// Step 4b - compute n_ind
	auto const n_ind = walk_size - last_enrichment();

	// Step 6 - Store the stats
	sW(indices) += W;
	Sn(indices) += 1;
	Se(indices) += (double) n_ind / walk_size;
      }
    }
  }

  void load(hdf5::handle const& loc) {
    std::cerr << "loading flatperm histograms: ";
    std::cerr << "sW, ";  hdf5::load(loc, sW, "sW");
    std::cerr << "Sn, ";  hdf5::load(loc, Sn, "Sn");
    std::cerr << "Se, ";  hdf5::load(loc, Se, "Se");
    std::cerr << "Enr, "; hdf5::load(loc, Enr, "Enr");
    std::cerr << "Pru\n"; hdf5::load(loc, Pru, "Pru");
  }

  void save(hdf5::handle const& loc) const {
    std::cerr << "saving flatperm histograms: ";
    std::cerr << "sW, ";  hdf5::save(loc, sW, "sW");
    std::cerr << "Sn, ";  hdf5::save(loc, Sn, "Sn");
    std::cerr << "Se, ";  hdf5::save(loc, Se, "Se");
    std::cerr << "Enr, "; hdf5::save(loc, Enr, "Enr");
    std::cerr << "Pru\n"; hdf5::save(loc, Pru, "Pru");
  }
};
}

#endif // FLATPERM_HPP

// vim: noai:ts=2:sw=2
