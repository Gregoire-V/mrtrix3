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

template <typename ValueType> class WeightedTransform : public TransformBase<ValueType> {
public:
  template <class MatrixType, class VectorType>
  WeightedTransform(const MatrixType &dirs, const VectorType &weights, int lmax, bool sym = true)
      : TransformBase<ValueType>(dirs, lmax, sym) {
    TransformBase<ValueType>::iSHT = wls(TransformBase<ValueType>::SHT, weights);
  }
};

} // namespace MR::Math::SH
