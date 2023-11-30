#pragma once
#include "./mode.hpp"

inline base_mode *create_mode(const std::string &mode_name) {
  if (mode_name == "full_factorization") return new ModeFullFactorization();
  if (mode_name == "full_factorization_pivots") return new ModeFullFactorizationPivots();
  if (mode_name == "vertex_factorization_pivots") return new ModeVertexFactorizationPivots();
  if (mode_name == "full_factorization_bath") return new ModeFullFactorizationBath();
  if (mode_name == "vertex_factorization") return new ModeVertexFactorization();
  if (mode_name == "tree_factorization1") return new ModeTreeFactorization1();
  if (mode_name == "tree_factorization1_bath") return new ModeTreeFactorization1Bath();
  if (mode_name == "vertex_factorization_symmetrized") return new ModeVertexFactorizationSymmetrized();
  if (mode_name == "vertex_factorization_bath") return new ModeVertexFactorizationBath();
  if (mode_name == "nested_tci") return new ModeNestedTCI();
  if (mode_name == "reuse_pivots") return new ModeReusePivots();
  if (mode_name == "partition_factorization") return new ModePartitionFactorization();
  if (mode_name == "combine_factorization") return new ModeCombineFactorization();
  if (mode_name == "explicit_sum") return new ModeExplicitSum();
  if (mode_name == "explicit_sum_exact") return new ModeExplicitSumExact();
  return nullptr;
}
