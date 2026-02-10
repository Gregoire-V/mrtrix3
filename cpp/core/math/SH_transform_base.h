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

template <typename ValueType> class TransformBase {
public:
  using matrix_type = Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic>;

  template <class MatrixType> TransformBase(const MatrixType &dirs, int lmax) : SHT(init_transform(dirs, lmax)) {}

  template <class VectorType1, class VectorType2> void A2SH(VectorType1 &sh, const VectorType2 &amplitudes) const {
    assert(iSHT.rows() > 0);
    sh.noalias() = iSHT * amplitudes;
  }
  template <class VectorType1, class VectorType2> void SH2A(VectorType1 &amplitudes, const VectorType2 &sh) const {
    amplitudes.noalias() = SHT * sh;
  }

  size_t n_SH() const { return SHT.cols(); }
  size_t n_amp() const { return SHT.rows(); }

  const matrix_type &mat_A2SH() const { return iSHT; }
  const matrix_type &mat_SH2A() const { return SHT; }

protected:
  matrix_type SHT, iSHT;
};

} // namespace MR::Math::SH
