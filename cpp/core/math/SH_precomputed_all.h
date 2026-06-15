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
#include "math/SH_precomputed_fraction.h"
#include "math/least_squares.h"
#include "math/legendre.h"
#include "mrtrix.h"

namespace MR::Math::SH {

// Needed to avoid an ambiguity between int and bool for PrecomputedAL constructors (otherwise doesn't compile)
struct SymFlag {
  SymFlag(bool sym) : value(sym) {}
  bool value;
};

//! Precomputed Associated Legrendre Polynomials - used to speed up SH calculation
template <typename ValueType> class PrecomputedAL {
public:
  using value_type = ValueType;

  // PrecomputedAL(bool sym = true) : lmax(0), ndir(0), nAL(0), inc(0.0), symmetric(sym) {}
  // PrecomputedAL(int up_to_lmax, bool sym) : symmetric(sym) { init(up_to_lmax, 512); }
  // PrecomputedAL(int up_to_lmax, int num_dir = 512, bool sym = true) : symmetric(sym) { init(up_to_lmax, num_dir); }

  PrecomputedAL(SymFlag sf = SymFlag(true)) : lmax(0), ndir(0), nAL(0), inc(0.0), symmetric(sf.value) {}
  PrecomputedAL(int up_to_lmax, SymFlag sf) : PrecomputedAL(up_to_lmax, 512, sf) {}
  PrecomputedAL(int up_to_lmax, int num_dir = 512, SymFlag sf = SymFlag(true)) : symmetric(sf.value) {
    init(up_to_lmax, num_dir);
  }

  bool operator!() const { return AL.empty(); }
  operator bool() const { return AL.size(); }

  void init(int up_to_lmax, int num_dir = 512) {
    lmax = up_to_lmax;
    ndir = num_dir;
    nAL = NforL_mpos(lmax, symmetric);
    inc = Math::pi / (ndir - 1);
    AL.resize(ndir * nAL);
    Eigen::Matrix<value_type, Eigen::Dynamic, 1, 0, 64> buf(lmax + 1);
    int l_step = (symmetric ? 2 : 1);

    for (int n = 0; n < ndir; n++) {
      typename std::vector<value_type>::iterator p = AL.begin() + n * nAL;
      value_type cos_el = std::cos(n * inc);
      for (int m = 0; m <= lmax; m++) {
        Legendre::Plm_sph(buf, lmax, m, cos_el);
        for (int l = (symmetric && (m & 1) ? m + 1 : m); l <= lmax; l += l_step)
          p[index_mpos(l, m, symmetric)] = buf[l];
      }
    }
  }

  void set(PrecomputedFraction<ValueType> &f, const ValueType inclination) const {
    f.f2 = inclination / inc;
    int i = static_cast<int>(std::trunc(f.f2));
    if (i < 0) {
      i = 0;
      f.f1 = 1.0;
      f.f2 = 0.0;
    } else if (i >= ndir - 1) {
      i = ndir - 1;
      f.f1 = 1.0;
      f.f2 = 0.0;
    } else {
      f.f2 -= i;
      f.f1 = 1.0 - f.f2;
    }

    f.p1 = AL.begin() + i * nAL;
    f.p2 = f.p1 + nAL;
  }

  ValueType get(const PrecomputedFraction<ValueType> &f, int i) const {
    ValueType v = f.f1 * f.p1[i];
    if (f.f2)
      v += f.f2 * f.p2[i];
    return v;
  }
  ValueType get(const PrecomputedFraction<ValueType> &f, int l, int m) const {
    return get(f, index_mpos(l, m, symmetric));
  }

  void get(ValueType *dest, const PrecomputedFraction<ValueType> &f) const {
    int l_step = (symmetric ? 2 : 1);
    for (int l = 0; l <= lmax; l += l_step) {
      for (int m = 0; m <= l; m++) {
        int i = index_mpos(l, m, symmetric);
        dest[i] = get(f, i);
      }
    }
  }

  template <class VectorType, class UnitVectorType>
  ValueType value(const VectorType &val, const UnitVectorType &unit_dir) const {
    PrecomputedFraction<ValueType> f;
    set(f, std::acos(unit_dir[2]));
    ValueType rxy = std::sqrt(pow2(unit_dir[1]) + pow2(unit_dir[0]));
    ValueType cp = (rxy) ? unit_dir[0] / rxy : 1.0;
    ValueType sp = (rxy) ? unit_dir[1] / rxy : 0.0;
    ValueType v = 0.0;
    int l_step = (symmetric ? 2 : 1);
    for (int l = 0; l <= lmax; l += l_step)
      v += get(f, l, 0) * val[index(l, 0, symmetric)];
    ValueType c0(1.0), s0(0.0);
    for (int m = 1; m <= lmax; m++) {
      ValueType c = c0 * cp - s0 * sp;
      ValueType s = s0 * cp + c0 * sp;
      for (int l = (symmetric && (m & 1) ? m + 1 : m); l <= lmax; l += l_step)
        v += get(f, l, m) * Math::sqrt2 * (c * val[index(l, m, symmetric)] + s * val[index(l, -m, symmetric)]);
      c0 = c;
      s0 = s;
    }
    return v;
  }

protected:
  int lmax, ndir, nAL;
  bool symmetric;
  ValueType inc;
  std::vector<ValueType> AL;
};

} // namespace MR::Math::SH
