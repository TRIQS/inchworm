#include <iostream>
#include <string>
#include "mode_factory.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  auto time_start = MPI_Wtime();

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
    NVTX_RANGE("main", 0);
    if (mode->rank == 0) { std::cout << "###### Tensor Train Based Hybridization Expansion Algorithm ######" << std::endl; }
    std::cout << std::fixed << std::setprecision(18);
    std::cerr << std::fixed << std::setprecision(18);
    mode->init(json_file_path, hyb_file_path);
    mode->run();
    if (mode->rank == 0) { mode->print_summary(); }
    delete mode;
  } else {
    std::cerr << "Invalid mode name" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  auto time = MPI_Wtime() - time_start;
  double time_sum;
  MPI_Reduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(mode->rank == 0) {
    std::cout << "CPU hours: " << time_sum / 3600 << std::endl;
  }
  MPI_Finalize();
  return 0;
}