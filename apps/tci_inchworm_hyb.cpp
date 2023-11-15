#include <iostream>
#include <string>
#include "mode/mode_factory.hpp"

int main(int argc, char *argv[]) {
  std::string mode_name {};
  std::string json_file_path {};

  // Check if mode and file path are provided
  if (argc > 2) {
    mode_name     = argv[1];
    json_file_path = argv[2];
  } else if (argc > 1) {
    mode_name     = argv[1];
    json_file_path = "../../apps/params.json";
  } else {
    std::cout << "Please provide mode (and file path)" << std::endl;
    std::cout << "Supported modes: full_factorization, vertex_factorization, nested_tci, reuse_pivots, partition_factorization, combine_factorization, explicit_sum" << std::endl;
    return 0;
  }

  base_mode *mode = create_mode(mode_name);
  if (mode) {
    std::cout << std::setprecision(18) << "##### Tensor Train Based Hybridyzation Expansion Inchworm Algorithm #####" << std::endl;
    mode->init(json_file_path);
    mode->prepare_input();
    mode->run_single_element();
    mode->print_summary();
    delete mode;
  } else {
    std::cout << "Invalid mode specified" << std::endl;
  }

  return 0;
}