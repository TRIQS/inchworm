// Copyright (c) 2026--present, The Simons Foundation
// This file is part of inchworm and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include <c2py/c2py.hpp>

#include <triqs/operators.hpp>
#include <triqs/atom_diag.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <nda/c2py/converters.hpp>
#include <triqs/c2py_converters/gf.hpp>
#include <triqs/c2py_converters/mesh.hpp>
#include <triqs/c2py_converters/fundamental_operator_set.hpp>
#include <triqs/c2py_converters/real_or_complex.hpp>
#include <triqs/c2py_converters/operators_real_complex.hpp>

#include <inchworm/atom_diag.hpp>

#include "ad_tools.wrap.cxx"
