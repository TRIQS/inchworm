#include <iostream>
#include <string>
#include "mode/mode_factory.hpp"

int main(int argc, char *argv[]) {
  std::string modeName;
  std::string jsonFilePath;

  // Check if mode and file path are provided
  if (argc > 2) {
    modeName         = argv[1];
    jsonFilePath = argv[2];
  } else if (argc > 1) {
    modeName         = argv[1];
    jsonFilePath = "../../apps/parameters.json";
  } else {
    std::cout << "Please provide mode (and file path)" << std::endl;
    return 0;
  }

  BaseMode *mode = createMode(modeName);
  if (mode) {
    mode->init(jsonFilePath);
    mode->prepareInput();
    mode->runSingleElement();
    delete mode; 
  } else {
    std::cout << "Invalid mode specified" << std::endl;
  }

  return 0;
}