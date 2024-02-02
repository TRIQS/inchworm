#include <iostream>
#include <string>
#include "mode/mode_factory.hpp"

int main(int argc, char *argv[]) {
  std::string mode_name {};
  std::string json_file_path {};

  // Check if mode and file path are provided
  if (argc == 3) {
    mode_name     = argv[1];
    json_file_path = argv[2];
  } else {
    std::cout << "Please provide mode name and parameter file" << std::endl;
    // TODO: Print usage
    return 0;
  }

  ModeBase *mode = create_mode(mode_name);
  if (mode) {
    std::cout << std::setprecision(18) << "##### Tensor Train Based Hybridization Expansion Algorithm #####" << std::endl;
    mode->init(json_file_path); // Read parameters from json file
    mode->run();
    mode->print_summary();
    delete mode;
  } else {
    std::cout << "Invalid mode specified" << std::endl;
  }

  return 0;
}