/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_KERNEL_LONG_GRAVITY_H
#define SWIFT_KERNEL_LONG_GRAVITY_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "approx_math.h"
#include "const.h"
#include "inline.h"

/* Standard headers */
#include <float.h>
#include <math.h>

#define GADGET2_LONG_RANGE_CORRECTION

#ifdef GADGET2_LONG_RANGE_CORRECTION
#define kernel_long_gravity_truncation_name "Gadget-like (using erfc())"
#else
#define kernel_long_gravity_truncation_name "Exp-based Sigmoid"
#endif

/**
 * @brief Derivatives of the long-range truncation function \f$\chi(r,r_s)\f$ up
 * to 5th order.
 */
struct chi_derivatives {

  /*! 0th order derivative \f$\chi(r,r_s)\f$ */
  double chi_0;

  /*! 1st order derivative \f$\partial_{r}\chi(r,r_s)\f$ */
  double chi_1;

  /*! 2nd order derivative \f$\partial_{rr}\chi(r,r_s)\f$ */
  double chi_2;

  /*! 3rd order derivative \f$\partial_{rrr}\chi(r,r_s)\f$ */
  double chi_3;

  /*! 4th order derivative \f$\partial_{rrrr}\chi(r,r_s)\f$ */
  double chi_4;

  /*! 5th order derivative \f$\partial_{rrrrr}\chi(r,r_s)\f$ */
  double chi_5;
};

/**
 * @brief Compute the derivatives of the long-range truncation function
 * \f$\chi(r,r_s)\f$ up to 5th order.
 *
 * @param r The distance.
 * @param r_s_inv The inverse of the long-range gravity mesh scale.
 * @param derivs (return) The computed #chi_derivatives.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_derivatives(
    const double r, const double r_s_inv, struct chi_derivatives *const derivs) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  /* Powers of u=r/2r_s */
  const double u = 0.5 * r * r_s_inv;
  const double u2 = u * u;
  const double u3 = u2 * u;
  const double u4 = u3 * u;

  /* Powers of (1/r_s) */
  const double r_s_inv2 = r_s_inv * r_s_inv;
  const double r_s_inv3 = r_s_inv2 * r_s_inv;
  const double r_s_inv4 = r_s_inv3 * r_s_inv;
  const double r_s_inv5 = r_s_inv4 * r_s_inv;

  /* Derivatives of \chi */
  derivs->chi_0 = approx_erfcf(u);
  derivs->chi_1 = -r_s_inv;
  derivs->chi_2 = r_s_inv2 * u;
  derivs->chi_3 = -r_s_inv3 * (u2 - 0.5);
  derivs->chi_4 = r_s_inv4 * (u3 - 1.5 * u);
  derivs->chi_5 = -r_s_inv5 * (u4 - 3. * u2 + 0.75);

  const double one_over_sqrt_pi = ((double)(M_2_SQRTPI * 0.5));
  const double common_factor = one_over_sqrt_pi * exp(-u2);

  /* Multiply in the common factors */
  derivs->chi_1 *= common_factor;
  derivs->chi_2 *= common_factor;
  derivs->chi_3 *= common_factor;
  derivs->chi_4 *= common_factor;
  derivs->chi_5 *= common_factor;

#else

  /* Powers of 2/r_s */
  const double c0 = 1;
  const double c1 = 2. * r_s_inv;
  const double c2 = c1 * c1;
  const double c3 = c2 * c1;
  const double c4 = c3 * c1;
  const double c5 = c4 * c1;

  /* 2r / r_s */
  const double x = c1 * r;

  /* e^(2r / r_s) */
  const double exp_x = exp(x);  // good_approx_expf(x);

  /* 1 / alpha(w) */
  const double a_inv = 1. + exp_x;

  /* Powers of alpha */
  const double a1 = 1. / a_inv;
  const double a2 = a1 * a1;
  const double a3 = a2 * a1;
  const double a4 = a3 * a1;
  const double a5 = a4 * a1;
  const double a6 = a5 * a1;

  /* Derivatives of \chi */
  derivs->chi_0 = -2. * exp_x * c0 * a1 + 2.;
  derivs->chi_1 = -2. * exp_x * c1 * a2;
  derivs->chi_2 = -2. * exp_x * c2 * (2. * a3 - a2);
  derivs->chi_3 = -2. * exp_x * c3 * (6. * a4 - 6. * a3 + a2);
  derivs->chi_4 = -2. * exp_x * c4 * (24. * a5 - 36. * a4 + 14. * a3 - a2);
  derivs->chi_5 = -2. * exp_x * c5 * (120. * a6 - 240. * a5 + 150. * a4 - 30. * a3 + a2);
#endif
}

/**
 * @brief Computes the long-range correction term for the potential calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_pot_eval(
    const float u, float *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float arg1 = u * 0.5f;
  const float term1 = approx_erfcf(arg1);

  *W = term1;
#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2 - 2 exp(x) * alpha */
  *W = 1.f - alpha * exp_x;
  *W *= 2.f;
#endif
}

/**
 * @brief Computes the long-range correction term for the force calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_force_eval(
    const float u, float *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));

  const float arg1 = u * 0.5f;
  const float arg2 = -arg1 * arg1;

  const float term1 = approx_erfcf(arg1);
  const float term2 = u * one_over_sqrt_pi * expf(arg2);

  *W = term1 + term2;
#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2*(x*alpha - x*alpha^2 - exp(x)*alpha + 1) */
  *W = 1.f - alpha;
  *W = *W * x - exp_x;
  *W = *W * alpha + 1.f;
  *W *= 2.f;
#endif
}

/**
 * @brief Returns the long-range truncation of the Poisson potential in Fourier
 * space.
 *
 * @param u2 The square of the Fourier mode times the cell scale
 * \f$u^2 = k^2r_s^2\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void fourier_kernel_long_grav_eval(
    const double u2, double *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION
  *W = exp(-u2);
#else
  const double u = sqrt(u2);
  const double arg = M_PI_2 * u;
  *W = arg / (sinh(arg) + FLT_MIN);
#endif
}

#endif /* SWIFT_KERNEL_LONG_GRAVITY_H */
