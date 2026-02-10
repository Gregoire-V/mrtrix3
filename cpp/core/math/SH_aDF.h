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

#include "math/SH_apodized_base.h"
#include "exception.h"
#include "math/least_squares.h"
#include "math/legendre.h"
#include "mrtrix.h"

namespace MR::Math::SH {

//! a class to hold the coefficients for an apodised disc function.
template <typename ValueType> class aDF : public ApodizedBase<ValueType> {
public:
  aDF(const size_t lmax) : ApodizedBase<ValueType>(lmax) {
    switch (lmax) {
    case 2:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.20980440;
      break;
    case 4:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.32978670, 0.07586349;
      break;
    case 6:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.39063368, 0.15768564, -0.03783370;
      break;
    case 8:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.42562849, 0.21421890, -0.08941999, 0.02256183;
      break;
    case 10:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.44556363, 0.25340653, -0.13204417, 0.05652718, -0.01499076;
      break;
    case 12:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.45894067, 0.28010314, -0.16643825, 0.08800988, -0.03868096,
          0.01075463;
      break;
    case 14:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.46756890, 0.29920102, -0.19231217, 0.11586324, -0.06198480,
          0.02794562, -0.00810281;
      break;
    case 16:
      ApodizedBase<ValueType>::RH << 1.00000000, -0.47395091, 0.31300714, -0.21252770, 0.13875942, -0.08428515,
          0.04562587, -0.02109019, 0.00635246;
      break;
    default:
      throw Exception("No aDF RH data for lmax " + str(lmax));
    }
  }
};

} // namespace MR::Math::SH
