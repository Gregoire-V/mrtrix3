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

//! a class to hold the coefficients for an apodised point-spread function.
template <typename ValueType> class aPSF : public ApodizedBase<ValueType> {
public:
  aPSF(const size_t lmax) : ApodizedBase<ValueType>(lmax) {
    switch (lmax) {
    case 2:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.41939279;
      break;
    case 4:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.63608543, 0.18487087;
      break;
    case 6:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.75490341, 0.37126442, 0.09614699;
      break;
    case 8:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.82384816, 0.51261696, 0.22440563, 0.05593079;
      break;
    case 10:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.86725945, 0.61519436, 0.34570667, 0.14300355, 0.03548062;
      break;
    case 12:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.89737759, 0.69278503, 0.45249879, 0.24169922, 0.09826171, 0.02502481;
      break;
    case 14:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.91717853, 0.74685644, 0.53467773, 0.33031863, 0.17013825, 0.06810155,
          0.01754930;
      break;
    case 16:
      ApodizedBase<ValueType>::RH << 1.00000000, 0.93261196, 0.79064858, 0.60562880, 0.41454703, 0.24880754, 0.12661242,
          0.05106681, 0.01365433;
      break;
    default:
      throw Exception("No aPSF RH data for lmax " + str(lmax));
    }
  }
};

} // namespace MR::Math::SH
