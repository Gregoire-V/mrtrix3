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

#include <Eigen/Dense>

#include "math/SH_aDF.h"
#include "math/SH_apodized_base.h"
#include "math/SH_aPSF.h"
#include "math/SH_functions.h"
#include "math/SH_helpers.h"
#include "math/SH_precomputed_all.h"
#include "math/SH_precomputed_fraction.h"
#include "math/SH_transform.h"
#include "math/SH_transform_base.h"
#include "math/SH_weighted_transform.h"
#include "math/least_squares.h"
#include "math/legendre.h"

namespace MR::Math::ZSH {

/** \defgroup zonal_spherical_harmonics Zonal Spherical Harmonics
 * \brief Classes & functions to manage zonal spherical harmonics
 * (spherical harmonic functions containing only m=0 terms). */

/** \addtogroup zonal_spherical_harmonics
 * @{ */

//! the number of (even-degree) coefficients for the given value of \a lmax
inline size_t NforL(int lmax) { return (1 + lmax / 2); }

//! compute the index for coefficient l
inline size_t index(int l) { return (l / 2); }

//! returns the largest \e lmax given \a N parameters
inline size_t LforN(int N) { return (2 * (N - 1)); }

} // namespace MR::Math::ZSH
