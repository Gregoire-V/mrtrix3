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

template <typename ValueType> class Transform {
public:
  using matrix_type = Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic>;

  template <class MatrixType>
  Transform(const MatrixType &dirs, const size_t lmax)
      // inclination angles are second column of azimuth/inclination matrix
      : ZSHT(init_amp_transform(dirs.col(1), lmax)), iZSHT(pinv(ZSHT)) {}

  template <class VectorType1, class VectorType2> void A2ZSH(VectorType1 &zsh, const VectorType2 &amplitudes) const {
    zsh.noalias() = iZSHT * amplitudes;
  }
  template <class VectorType1, class VectorType2> void ZSH2A(VectorType1 &amplitudes, const VectorType2 &zsh) const {
    amplitudes.noalias() = ZSHT * zsh;
  }

  size_t n_ZSH() const { return ZSHT.cols(); }
  size_t n_amp() const { return ZSHT.rows(); }

  const matrix_type &mat_A2ZSH() const { return iZSHT; }
  const matrix_type &mat_ZSH2A() const { return ZSHT; }

protected:
  matrix_type ZSHT, iZSHT;
}; 

} // namespace MR::Math::ZSH
