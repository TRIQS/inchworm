#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "./mode.hpp"

using namespace inchworm;

void ModeDebug::run() {
  // debug model now only support using exact propagator to evalute a single element of the propagator
  if (gp.target != "propagator" || gp.model_type != 0) {
    std::cerr << "debug mode: not implemented yet" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  
}
