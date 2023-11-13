#include "./mode.hpp"

inline BaseMode* createMode(const std::string& modeName) {
    if (modeName == "use_norm_pivots") return new ModeUseNormPivots();
    return nullptr;
}
