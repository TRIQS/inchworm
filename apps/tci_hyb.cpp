#include <iostream>
#include <string>
#include "mode/mode_factory.hpp"

int main(int argc, char *argv[]) {
  std::string mode_name{};
  std::string json_file_path{};
  std::string hyb_file_path{};

  // Check if mode and file path are provided
  if (argc == 3) {
    mode_name      = argv[1];
    json_file_path = argv[2];
  } 
  else if (argc == 4) {
    mode_name      = argv[1];
    json_file_path = argv[2];
    hyb_file_path= argv[3];
  }
  else {
    std::cout << "Please provide mode name and parameter file" << std::endl;
    std::cout << "Usage: tci_hyb <mode> <json_file_path>" << std::endl;
    std::cout << "mode: inch, debug, bare" << std::endl;
    std::cout << "--- global parameters ---" << std::endl;
    std::cout << "target: propagator, green_function, both" << std::endl;
    std::cout << "integrand: plain, sum_phi" << std::endl;
    std::cout << "integral_variable: v, v_iota, v_id, v_iota_id" << std::endl;
    std::cout << "tci_shape: plain, partition, vertex" << std::endl;
    std::cout << "trick: none, scan_pivots, pretraining" << std::endl;
    std::cout << "model_type: 0 (discrete bath), 1 (read hybridization from file), 2(Bethe lattice)" << std::endl;
  }

  ModeBase *mode = create_mode(mode_name);
  if (mode) {
    std::cout << std::setprecision(18) << "###### Tensor Train Based Hybridization Expansion Algorithm ######" << std::endl;
    mode->init(json_file_path,hyb_file_path);
    mode->run();
    mode->print_summary();
    delete mode;
  } else {
    std::cout << "Invalid mode specified" << std::endl;
  }

  return 0;
}