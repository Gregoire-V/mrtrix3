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
#include "math/SH_helpers.h"
#include "math/SH_precomputed_all.h"
#include "math/least_squares.h"
#include "math/legendre.h"
#include "mrtrix.h"

namespace MR::Math::SH {

//! form the SH->amplitudes matrix
/*! This computes the matrix \a SHT mapping spherical harmonic
 * coefficients up to maximum harmonic degree \a lmax onto directions \a
 * dirs (in spherical coordinates, with columns [ azimuth inclination ]). */
template <class MatrixType>
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>
init_transform(const MatrixType &dirs, const int lmax, bool sym = true) {
  using namespace Eigen;
  using value_type = typename MatrixType::Scalar;
  ssize_t l_step = (sym ? 2 : 1);
  if (dirs.cols() != 2)
    throw Exception("direction matrix should have 2 columns: [ azimuth inclination ]");
  Matrix<value_type, Dynamic, Dynamic> SHT(dirs.rows(), NforL(lmax, sym));
  Matrix<value_type, Dynamic, 1, 0, 64> AL(lmax + 1);
  for (ssize_t i = 0; i < dirs.rows(); i++) {
    const value_type z = std::cos(dirs(i, 1));
    Legendre::Plm_sph(AL, lmax, 0, z);
    for (int l = 0; l <= lmax; l += l_step)
      SHT(i, index(l, 0, sym)) = AL[l];
    for (int m = 1; m <= lmax; m++) {
      Legendre::Plm_sph(AL, lmax, m, z);
      for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step) {
        SHT(i, index(l, m, sym)) = Math::sqrt2 * AL[l] * std::cos(m * dirs(i, 0));
        SHT(i, index(l, -m, sym)) = Math::sqrt2 * AL[l] * std::sin(m * dirs(i, 0));
      }
    }
  }
  return SHT;
}

//! form the SH->amplitudes matrix
/*! This computes the matrix \a SHT mapping spherical harmonic
 * coefficients up to maximum harmonic degree \a lmax onto directions \a
 * dirs (in cartesian coordinates, with columns [ x y z ] ans normalised). */
template <class MatrixType>
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>
init_transform_cart(const MatrixType &dirs, const int lmax, bool sym = true) {
  using namespace Eigen;
  using value_type = typename MatrixType::Scalar;
  ssize_t l_step = (sym ? 2 : 1);
  if (dirs.cols() != 3)
    throw Exception("direction matrix should have 3 columns: [ x y z ]");
  Matrix<value_type, Dynamic, Dynamic> SHT(dirs.rows(), NforL(lmax, sym));
  Matrix<value_type, Dynamic, 1, 0, 64> AL(lmax + 1);
  for (ssize_t i = 0; i < dirs.rows(); i++) {
    value_type z = dirs(i, 2);
    value_type rxy = std::hypot(dirs(i, 0), dirs(i, 1));
    value_type cp = (rxy) ? dirs(i, 0) / rxy : 1.0;
    value_type sp = (rxy) ? dirs(i, 1) / rxy : 0.0;
    Legendre::Plm_sph(AL, lmax, 0, z);
    for (int l = 0; l <= lmax; l += l_step)
      SHT(i, index(l, 0, sym)) = AL[l];
    value_type c0(1.0), s0(0.0);
    for (int m = 1; m <= lmax; m++) {
      Legendre::Plm_sph(AL, lmax, m, z);
      value_type c = c0 * cp - s0 * sp;
      value_type s = s0 * cp + c0 * sp;
      for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step) {
        SHT(i, index(l, m, sym)) = Math::sqrt2 * AL[l] * c;
        SHT(i, index(l, -m, sym)) = Math::sqrt2 * AL[l] * s;
      }
      c0 = c;
      s0 = s;
    }
  }
  return SHT;
}

//! scale the coefficients of each SH degree by the corresponding value in \a coefs
template <class MatrixType, class VectorType>
inline void scale_degrees_forward(MatrixType &SH2amp_mapping, const VectorType &coefs, bool sym = true) {
  ssize_t l = 0, nl = 1;
  ssize_t l_factor = (sym ? 2 : 1);
  for (ssize_t col = 0; col < SH2amp_mapping.cols(); ++col) {
    if (col >= nl) {
      l++;
      nl = NforL(l_factor * l, sym);
    }
    SH2amp_mapping.col(col) *= coefs[l];
  }
}

//! scale the coefficients of each SH degree by the corresponding value in \a coefs
template <typename MatrixType, class VectorType>
inline void scale_degrees_inverse(MatrixType &amp2SH_mapping, const VectorType &coefs, bool sym = true) {
  ssize_t l = 0, nl = 1;
  ssize_t l_factor = (sym ? 2 : 1);
  for (ssize_t row = 0; row < amp2SH_mapping.rows(); ++row) {
    if (row >= nl) {
      l++;
      nl = NforL(l_factor * l, sym);
    }
    amp2SH_mapping.row(row) *= coefs[l];
  }
}

//! invert any non-zero coefficients in \a coefs
template <typename VectorType>
inline Eigen::Matrix<typename VectorType::Scalar, Eigen::Dynamic, 1> invert(const VectorType &coefs) {
  Eigen::Matrix<typename VectorType::Scalar, Eigen::Dynamic, 1> ret(coefs.size());
  for (size_t n = 0; n < coefs.size(); ++n)
    ret[n] = (coefs[n] ? 1.0 / coefs[n] : 0.0);
  return ret;
}

template <class VectorType>
inline typename VectorType::Scalar value(const VectorType &coefs,
                                         typename VectorType::Scalar cos_inclination,
                                         typename VectorType::Scalar cos_azimuth,
                                         typename VectorType::Scalar sin_azimuth,
                                         int lmax,
                                         bool sym = true) {
  using value_type = typename VectorType::Scalar;
  ssize_t l_step = (sym ? 2 : 1);
  value_type amplitude = 0.0;
  Eigen::Matrix<value_type, Eigen::Dynamic, 1, 0, 64> AL(lmax + 1);
  Legendre::Plm_sph(AL, lmax, 0, cos_inclination);
  for (int l = 0; l <= lmax; l += l_step)
    amplitude += AL[l] * coefs[index(l, 0, sym)];
  value_type c0(1.0), s0(0.0);
  for (int m = 1; m <= lmax; m++) {
    Legendre::Plm_sph(AL, lmax, m, cos_inclination);
    value_type c = c0 * cos_azimuth - s0 * sin_azimuth; // std::cos(m*azimuth)
    value_type s = s0 * cos_azimuth + c0 * sin_azimuth; // std::sin(m*azimuth)
    for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step)
      amplitude += AL[l] * Math::sqrt2 * (c * coefs[index(l, m, sym)] + s * coefs[index(l, -m, sym)]);
    c0 = c;
    s0 = s;
  }
  return amplitude;
}

template <class VectorType>
inline typename VectorType::Scalar value(const VectorType &coefs,
                                         typename VectorType::Scalar cos_inclination,
                                         typename VectorType::Scalar azimuth,
                                         int lmax,
                                         bool sym = true) {
  return value(coefs, cos_inclination, std::cos(azimuth), std::sin(azimuth), lmax, sym);
}

template <class VectorType1, class VectorType2>
inline typename VectorType1::Scalar
value(const VectorType1 &coefs, const VectorType2 &unit_dir, int lmax, bool sym = true) {
  using value_type = typename VectorType1::Scalar;
  value_type rxy = std::sqrt(pow2(unit_dir[1]) + pow2(unit_dir[0]));
  value_type cp = (rxy) ? unit_dir[0] / rxy : 1.0;
  value_type sp = (rxy) ? unit_dir[1] / rxy : 0.0;
  return value(coefs, unit_dir[2], cp, sp, lmax, sym);
}

template <class VectorType1, class VectorType2>
inline VectorType1 &delta(VectorType1 &delta_vec, const VectorType2 &unit_dir, int lmax, bool sym = true) {
  using value_type = typename VectorType1::Scalar;
  ssize_t l_step = (sym ? 2 : 1);
  delta_vec.resize(NforL(lmax, sym));
  value_type rxy = std::sqrt(pow2(unit_dir[1]) + pow2(unit_dir[0]));
  value_type cp = (rxy) ? unit_dir[0] / rxy : 1.0;
  value_type sp = (rxy) ? unit_dir[1] / rxy : 0.0;
  Eigen::Matrix<value_type, Eigen::Dynamic, 1, 0, 64> AL(lmax + 1);
  Legendre::Plm_sph(AL, lmax, 0, unit_dir[2]);
  for (int l = 0; l <= lmax; l += l_step)
    delta_vec[index(l, 0, sym)] = AL[l];
  value_type c0(1.0), s0(0.0);
  for (int m = 1; m <= lmax; m++) {
    Legendre::Plm_sph(AL, lmax, m, unit_dir[2]);
    value_type c = c0 * cp - s0 * sp;
    value_type s = s0 * cp + c0 * sp;
    for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step) {
      delta_vec[index(l, m, sym)] = AL[l] * Math::sqrt2 * c;
      delta_vec[index(l, -m, sym)] = AL[l] * Math::sqrt2 * s;
    }
    c0 = c;
    s0 = s;
  }
  return delta_vec;
}

template <class VectorType1, class VectorType2> inline VectorType1 &SH2RH(VectorType1 &RH, const VectorType2 &sh) {
  using value_type = typename VectorType2::Scalar;
  RH.resize(sh.size());
  int lmax = 2 * sh.size() + 1;
  Eigen::Matrix<value_type, Eigen::Dynamic, 1, 0, 64> AL(lmax + 1);
  Legendre::Plm_sph(AL, lmax, 0, 1.0);
  for (ssize_t l = 0; l < sh.size(); l++)
    RH[l] = sh[l] / AL[2 * l];
  return RH;
}

template <class VectorType>
inline Eigen::Matrix<typename VectorType::Scalar, Eigen::Dynamic, 1> SH2RH(const VectorType &sh) {
  Eigen::Matrix<typename VectorType::Scalar, Eigen::Dynamic, 1> RH(sh.size());
  SH2RH(RH, sh);
  return RH;
}

//! perform spherical convolution, in place
/*! perform spherical convolution of SH coefficients \a sh with response
 * function \a RH, storing the results in place in vector \a sh. */
template <class VectorType1, class VectorType2> inline VectorType1 &sconv(VectorType1 &sh, const VectorType2 &RH) {
  assert(static_cast<size_t>(sh.size()) >= NforL(2 * (RH.size() - 1)));
  for (ssize_t i = 0; i < RH.size(); ++i) {
    int l = 2 * i;
    for (int m = -l; m <= l; ++m)
      sh[index(l, m)] *= RH[i];
  }
  return sh;
}

//! perform spherical convolution
/*! perform spherical convolution of SH coefficients \a sh with response
 * function \a RH, storing the results in vector \a C. */
template <class VectorType1, class VectorType2, class VectorType3>
inline VectorType1 &sconv(VectorType1 &C, const VectorType2 &RH, const VectorType3 &sh) {
  assert(static_cast<size_t>(sh.size()) >= NforL(2 * (RH.size() - 1)));
  C.resize(NforL(2 * (RH.size() - 1)));
  for (ssize_t i = 0; i < RH.size(); ++i) {
    int l = 2 * i;
    for (int m = -l; m <= l; ++m)
      C[index(l, m)] = RH[i] * sh[index(l, m)];
  }
  return C;
}

//! perform spherical convolution, in place
/*! perform spherical convolution of SH coefficients, stored in rows
 * in matrix \a sh with response function \a RH, storing the results
 * in place in matrix \a sh. */
template <class MatrixType1, class VectorType2> inline MatrixType1 &sconv_mat(MatrixType1 &sh, const VectorType2 &RH) {
  assert(static_cast<size_t>(sh.cols()) >= NforL(2 * (RH.size() - 1)));
  for (ssize_t i = 0; i < RH.size(); ++i) {
    int l = 2 * i;
    for (int m = -l; m <= l; ++m)
      sh.col(index(l, m)) *= RH[i];
  }
  return sh;
}

//! estimate direction & amplitude of SH peak
/*! find a peak of an SH series using Gauss-Newton optimisation, modified
 * to operate directly in spherical coordinates. The initial search
 * direction is \a unit_init_dir. If \a precomputer is not nullptr, it
 * will be used to speed up the calculations, at the cost of a minor
 * reduction in accuracy. */
template <class VectorType, class UnitVectorType>
inline typename VectorType::Scalar get_peak(const VectorType &sh,
                                            int lmax,
                                            UnitVectorType &unit_init_dir,
                                            PrecomputedAL<typename VectorType::Scalar> *precomputer = nullptr,
                                            bool sym = true) {
  static const default_type max_dir_change = 0.2;
  static const default_type angle_tolerance = 1e-4;
  using value_type = typename VectorType::Scalar;
  assert(std::isfinite(unit_init_dir[0]));
  for (int i = 0; i < 50; i++) {
    value_type az = std::atan2(unit_init_dir[1], unit_init_dir[0]);
    value_type el = std::acos(unit_init_dir[2]);
    value_type amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2;
    derivatives(sh, lmax, el, az, amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2, precomputer, sym);

    value_type del = sqrt(dSH_del * dSH_del + dSH_daz * dSH_daz);
    value_type daz = 0.0;
    if (del != 0.0) {
      daz = dSH_daz / del;
      del = dSH_del / del;
    }

    value_type dSH_dt = daz * dSH_daz + del * dSH_del;
    value_type d2SH_dt2 = daz * daz * d2SH_daz2 + 2.0 * daz * del * d2SH_deldaz + del * del * d2SH_del2;
    value_type dt = d2SH_dt2 ? (-dSH_dt / d2SH_dt2) : 0.0;

    if (dt < 0.0)
      dt = -dt;
    if (dt > max_dir_change)
      dt = max_dir_change;

    del *= dt;
    daz *= dt;

    unit_init_dir[0] += del * std::cos(az) * std::cos(el) - daz * std::sin(az);
    unit_init_dir[1] += del * std::sin(az) * std::cos(el) + daz * std::cos(az);
    unit_init_dir[2] -= del * std::sin(el);
    unit_init_dir.normalize();

    if (dt < angle_tolerance)
      return amplitude;
  }

  unit_init_dir.fill(std::numeric_limits<typename UnitVectorType::Scalar>::quiet_NaN());
  DEBUG("failed to find SH peak!");
  return std::numeric_limits<typename VectorType::Scalar>::quiet_NaN();
}

//! computes first and second order derivatives of SH series
/*! This is used primarily in the get_peak() function. */
template <class VectorType>
inline void derivatives(const VectorType &sh,
                        const int lmax,
                        const typename VectorType::Scalar inclination,
                        const typename VectorType::Scalar azimuth,
                        typename VectorType::Scalar &amplitude,
                        typename VectorType::Scalar &dSH_del,
                        typename VectorType::Scalar &dSH_daz,
                        typename VectorType::Scalar &d2SH_del2,
                        typename VectorType::Scalar &d2SH_deldaz,
                        typename VectorType::Scalar &d2SH_daz2,
                        PrecomputedAL<typename VectorType::Scalar> *precomputer,
                        bool sym = true) {
  if (lmax < 0) {
    throw std::logic_error("lmax cannot be negative!");
  }

  using value_type = typename VectorType::Scalar;
  const value_type sin_incl = std::sin(inclination);
  const value_type cos_incl = std::cos(inclination);
  bool atpole = sin_incl < 1e-4;
  ssize_t l_step = (sym ? 2 : 1);

  dSH_del = dSH_daz = d2SH_del2 = d2SH_deldaz = d2SH_daz2 = 0.0;
  VLA_MAX(AL, value_type, NforL_mpos(lmax, sym), 64);

  if (precomputer) {
    PrecomputedFraction<value_type> f;
    precomputer->set(f, inclination);
    precomputer->get(AL, f);
  } else {
    Eigen::Matrix<value_type, Eigen::Dynamic, 1, 0, 64> buf(lmax + 1);
    for (int m = 0; m <= lmax; m++) {
      Legendre::Plm_sph(buf, lmax, m, cos_incl);
      for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step)
        AL[index_mpos(l, m, sym)] = buf[l];
    }
  }

  amplitude = sh[index(0, 0, sym)] * AL[index_mpos(0, 0, sym)];
  for (int l = l_step; l <= lmax; l += l_step) {
    const value_type &v(sh[index(l, 0, sym)]);
    amplitude += v * AL[index_mpos(l, 0, sym)];
    dSH_del += v * sqrt(static_cast<value_type>(l * (l + 1))) * AL[index_mpos(l, 1, sym)];
    d2SH_del2 += v *
                 (sqrt(static_cast<value_type>(l * (l + 1) * (l - 1) * (l + 2))) * AL[index_mpos(l, 2, sym)] -
                  l * (l + 1) * AL[index_mpos(l, 0, sym)]) /
                 2.0;
  }

  for (int m = 1; m <= lmax; m++) {
    value_type caz = Math::sqrt2 * std::cos(m * azimuth);
    value_type saz = Math::sqrt2 * std::sin(m * azimuth);
    for (int l = (sym && (m & 1) ? m + 1 : m); l <= lmax; l += l_step) {
      const value_type &vp(sh[index(l, m, sym)]);
      const value_type &vm(sh[index(l, -m, sym)]);
      amplitude += (vp * caz + vm * saz) * AL[index_mpos(l, m, sym)];

      value_type tmp = sqrt(static_cast<value_type>((l + m) * (l - m + 1))) * AL[index_mpos(l, m - 1, sym)];
      if (l > m)
        tmp -= sqrt(static_cast<value_type>((l - m) * (l + m + 1))) * AL[index_mpos(l, m + 1, sym)];
      tmp /= -2.0;
      dSH_del += (vp * caz + vm * saz) * tmp;

      value_type tmp2 = -((l + m) * (l - m + 1) + (l - m) * (l + m + 1)) * AL[index_mpos(l, m, sym)];
      if (m == 1)
        tmp2 -= sqrt(static_cast<value_type>((l + m) * (l - m + 1) * (l + m - 1) * (l - m + 2))) *
                AL[index_mpos(l, 1, sym)];
      else
        tmp2 += sqrt(static_cast<value_type>((l + m) * (l - m + 1) * (l + m - 1) * (l - m + 2))) *
                AL[index_mpos(l, m - 2, sym)];
      if (l > m + 1)
        tmp2 += sqrt(static_cast<value_type>((l - m) * (l + m + 1) * (l - m - 1) * (l + m + 2))) *
                AL[index_mpos(l, m + 2, sym)];
      tmp2 /= 4.0;
      d2SH_del2 += (vp * caz + vm * saz) * tmp2;

      if (atpole)
        dSH_daz += (vm * caz - vp * saz) * tmp;
      else {
        d2SH_deldaz += m * (vm * caz - vp * saz) * tmp;
        dSH_daz += m * (vm * caz - vp * saz) * AL[index_mpos(l, m, sym)];
        d2SH_daz2 -= (vp * caz + vm * saz) * m * m * AL[index_mpos(l, m, sym)];
      }
    }
  }

  if (!atpole) {
    dSH_daz /= sin_incl;
    d2SH_deldaz /= sin_incl;
    d2SH_daz2 /= sin_incl * sin_incl;
  }
}

//! determines whether the input is symmetric or asymmetric
/*! Let N be the number of SH coefficients. For some values of N, there exists both
 * l_sym such that N = NforL(l_sym, sym = True) and l_asym such that
 * N = NforL(l_asym, sym = False). The first two are N = 1 and N = 36. If N = 1, we
 * know that the input is symmetric. If N = 36, l_asym = 5 and l_sym = 7. But as
 * l_sym is odd while it should be even, we consider the input to be asymmetric.
 * The next ambiguity appears when N = 1225 (l_asym = 34, l_sym = 48). Hence, the
 * output of this function can be considered valid for N < 1225, or always valid if
 * we assume that such a high value will never be encountered. */
template <class ImageType> bool is_symetric(const ImageType &H) {
  bool symmetric_input = true;
  size_t n_SH = H.size(3);
  double sqrt_n_SH = std::sqrt(n_SH);
  if (n_SH != 1) {
    // the input is asymmetric if and only if n_SH is a perfect square
    symmetric_input = static_cast<size_t>(std::floor(sqrt_n_SH)) != static_cast<size_t>(std::ceil(sqrt_n_SH));
  }
  return symmetric_input;
}

//! convenience function to check if an input image can contain SH coefficients
template <class ImageType> void check(const ImageType &H, bool authorize_asym = false) {
  if (H.ndim() < 4)
    throw Exception("image \"" + H.name() + "\" does not contain SH coefficients - not 4D");
  if (!authorize_asym) {
    size_t l = LforN(H.size(3));
    if (l % 2 || NforL(l) != static_cast<size_t>(H.size(3)))
      throw Exception("image \"" + H.name() +
                      "\" does not contain SH coefficients - unexpected number of coefficients");
  } else {
    if (!feasible_N(H.size(3), true) && !feasible_N(H.size(3), false))
      throw Exception("image \"" + H.name() +
                      "\" does not contain SH coefficients - unexpected number of coefficients");
  }
}
/** @} */

} // namespace MR::Math::SH
