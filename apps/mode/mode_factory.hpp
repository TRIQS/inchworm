#include "./mode.hpp"

inline BaseMode* createMode(const std::string& modeName) {
    if (modeName == "use_norm_pivots") return new ModeUseNormPivots();
    if (modeName == "full_factorization") return new ModeFullFactorization();
    if (modeName == "vertex_factorization") return new ModeVertexFactorization();
    if (modeName == "nested_tci") return new ModeNestedTCI();
    if (modeName == "reuse_pivots") return new ModeReusePivots();
    if (modeName == "partition_factorization") return new ModePartitionFactorization();
    if (modeName == "combine_factorization") return new ModeCombineFactorization();
    if (modeName == "explicit_sum") return new ModeExplicitSum();
    return nullptr;
}
