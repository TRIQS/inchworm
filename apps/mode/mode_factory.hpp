#pragma once
#include "./mode.hpp"

inline ModeBase *create_mode(const std::string &mode_name) {
  if (mode_name == "debug") return new ModeDebug();
  if (mode_name == "inch") return new ModeInchworm();
  if (mode_name == "bare") return new ModeBare();
  return nullptr;
}
