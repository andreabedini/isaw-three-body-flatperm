/*
 * main.cpp
 *
 */

#include "instance.hpp"

#include "hdf5pp/hdf5.hpp"

#include "hdf5_hl.h"

#include <boost/asio/signal_set.hpp>
#include <boost/asio/deadline_timer.hpp>
#include <boost/cstdint.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>

#include <csignal>
#include <string>
#include <stdexcept>

using boost::uint64_t;

struct handlers {
  std::function<void()> callback;

  boost::asio::signal_set signals;
  boost::asio::deadline_timer timer;
  boost::posix_time::time_duration const checkpoint_interval;

  handlers(boost::asio::io_service& io_service, std::function<void()> callback)
    : callback(callback)
    , signals(io_service, SIGHUP, SIGTERM, SIGINT)
    , timer(io_service)
    , checkpoint_interval(boost::posix_time::hours(1))
  {
    signals.async_wait(boost::bind(&handlers::do_signal, this, _1, _2));
    // set the timer
    timer.expires_from_now(checkpoint_interval);
    timer.async_wait(boost::bind(&handlers::do_timer, this, _1));
  }

  void do_signal(const boost::system::error_code&,
                 int signal_number)
  {
    if (signal_number == SIGHUP) {
      callback();
      signals.async_wait(boost::bind(&handlers::do_signal, this, _1, _2));
    } else
      signals.get_io_service().stop();
  }

  void do_timer(const boost::system::error_code&)
  {
    callback();

    // reset the timer
    timer.expires_from_now(checkpoint_interval);
    timer.async_wait(boost::bind(&handlers::do_timer, this, _1));
  }
};

int main(int argc, char* argv[])
{
  std::cerr << "this is " << PACKAGE << "\n";

  namespace po = boost::program_options;
  namespace pt = boost::posix_time;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")

    ("resume", "resume")

    ("filename",        po::value<std::string>()->required(),
     "filename (required)")

    ("tours,S",         po::value<unsigned int>(),
     "tours")

    ("seed",            po::value<unsigned int>(),
     "random generator seed")

    ("length,N",        po::value<unsigned int>(),
     "walk maximum length (ignored when resuming)")

    ("mu",              po::value<double>()->default_value(1),
     "weight renormalization  (ignored when resuming)")
    ;

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (po::required_option e) {
    std::cerr << desc << "\n";
    return 0;
  }

  if (vm.count("help")) {
    std::cerr << desc << "\n";
    return 0;
  }

  const std::string filename = vm["filename"].as<std::string>();
  const unsigned int S = vm["tours"].as<unsigned int>();

  //////////////////////////////////////////////////

  hdf5::file hfile = vm.count("resume")
    ? hdf5::file::open  (filename, H5F_ACC_RDWR)
    : hdf5::file::create(filename, H5F_ACC_TRUNC)
    ;

  instance my_instance = vm.count("resume")
    ? instance(hfile)
    : instance(vm["length"].as<unsigned int>(),
               vm["mu"].as<double>())
    ;

  //
  // default seed is 1
  //
  unsigned int seed = vm.count("seed")
    ? vm["seed"].as<unsigned int>()
    : 1
    ;

  std::cerr << "seed set to " << seed << "\n";
  my_instance.rng.seed(seed);

  //////////////////////////////////////////////////////////////////////

  boost::asio::io_service io_service;

  auto save_data = [&] {
    my_instance.save(hfile);
    hfile.flush();
    my_instance.print_stats();
  };

  handlers my_handlers{io_service, save_data};

  //////////////////////////////////////////////////
  boost::thread t([&] {
    try {
      my_instance.run(S);
      io_service.stop();
    } catch (boost::thread_interrupted e) {
      std::cerr << "interrupted!\n";
    }
    });

  io_service.run();
  std::cerr << "main thread ready to stop.\n";
  t.interrupt();
  t.join();
  save_data();
}
