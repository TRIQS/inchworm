#pragma once
#include "./mode.hpp"

inline base_mode* create_mode(const std::string& mode_name) {
    if (mode_name == "full_factorization") return new ModeFullFactorization();
    if (mode_name == "vertex_factorization") return new ModeVertexFactorization();
    if (mode_name == "nested_tci") return new ModeNestedTCI();
    if (mode_name == "reuse_pivots") return new ModeReusePivots();
    if (mode_name == "partition_factorization") return new ModePartitionFactorization();
    if (mode_name == "combine_factorization") return new ModeCombineFactorization();
    if (mode_name == "explicit_sum") return new ModeExplicitSum();
    return nullptr;
}
