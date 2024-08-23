#include <iostream>
#include <string>
#include "mode_factory.hpp"

int main(int argc, char *argv[]) {
  std::string mode_name{};
  std::string json_file_path{};
  std::string hyb_file_path{};
  // Check if mode and file path are provided
  if (argc == 3) {
    mode_name      = argv[1];
    json_file_path = argv[2];
  } else if (argc == 4) {
    mode_name      = argv[1];
    json_file_path = argv[2];
    hyb_file_path  = argv[3];
  } else {
    std::cout << "Please provide mode name and parameter file" << std::endl;
    std::cout << "Usage: tci_hyb <mode> <json_file_path>" << std::endl;
    std::cout << "mode: inch, debug, bare" << std::endl;
    std::cout << "--- parameters in the json file ---" << std::endl;
  }

  ModeBase *mode = create_mode(mode_name);
  if (mode) {
    std::cout << "###### Tensor Train Based Hybridization Expansion Algorithm ######" << std::endl;
    std::cerr << std::setprecision(18);
    mode->init(json_file_path, hyb_file_path);
    double start_time = omp_get_wtime();
    mode->run();
    std::cout << "time: ";
    double end_time     = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    std::cout << elapsed_time << " sec" << std::endl;
    mode->print_summary();
    delete mode;
  } else {
    std::cout << "Invalid mode specified" << std::endl;
  }

  return 0;
}