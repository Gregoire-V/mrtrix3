/* Copyright (c) 2008-2026 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Covered Software is provided under this License on an "as is"
 * basis, without warranty of any kind, either expressed, implied, or
 * statutory, including, without limitation, warranties that the
 * Covered Software is free of defects, merchantable, fit for a
 * particular purpose or non-infringing.
 * See the Mozilla Public License v. 2.0 for more details.
 *
 * For more details, see http://www.mrtrix.org/.
 */

#pragma once

#include <string>

#include "exception.h"
#include "math/least_squares.h"
#include "math/legendre.h"
#include "mrtrix.h"

namespace MR::Math::SH {

/** \defgroup spherical_harmonics Spherical Harmonics
 * \brief Classes & functions to manage spherical harmonics. */

/** \addtogroup spherical_harmonics
 * @{ */

//! a string containing a description of the SH storage convention
/*! This can used directly in the DESCRIPTION field of a command's
 * usage() function. */
extern const std::string encoding_description;

//! the number of coefficients for the given value of \a lmax
inline size_t NforL(int lmax, bool sym = true) { return sym ? (lmax + 1) * (lmax + 2) / 2 : (lmax + 1) * (lmax + 1); }

//! compute the index for coefficient (l,m)
inline size_t index(int l, int m, bool sym = true) { return sym ? l * (l + 1) / 2 + m : l * (l + 1) + m; }

//! same as NforL(), but consider only non-negative orders \e m
inline size_t NforL_mpos(int lmax, bool sym = true) { return sym ? (lmax / 2 + 1) * (lmax / 2 + 1) : (lmax + 1) * (lmax + 2) / 2; }

//! same as index(), but consider only non-negative orders \e m
inline size_t index_mpos(int l, int m, bool sym = true) { return sym ? l * l / 4 + m : l * (l + 1) / 2 + m; }

//! returns the largest \e lmax given \a N parameters
inline size_t LforN(int N, bool sym = true) { 
  if (sym)
  {
    return N ? 2 * std::floor<size_t>((std::sqrt(static_cast<default_type>(1 + 8 * N)) - 3.0) / 4.0) : 0;
  }
  else
  {
    return N ? std::floor(std::sqrt(N)) - 1 : 0;
  } 
}

//! returns whether a cardinality is commensurate with a set of SH coefficients
inline bool feasible_N(int N, bool sym = true) { return NforL(LforN(N, sym), sym) == N; }

/* Never used, so could be removed
namespace {
template <typename> struct __dummy {
  using type = int;
};
} // namespace
 */

/** @} */
} // namespace MR::Math::SH
