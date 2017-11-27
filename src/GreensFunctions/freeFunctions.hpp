#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "DefsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

double expxsq_erfc(double x);

double W(double a, double b);

/* Functions needed in 1D Green's functions. For explanation, see the cpp file. */

double XP00(double r, double t, double r0, double D, double v);

double XI00(double r, double t, double r0, double D, double v);

//double XS00(double t, double r0, double D, double v);

double XP10(double r, double t, double r0, double D, double v);

double XI10(double r, double t, double r0, double D, double v);

double XS10(double t, double r0, double D, double v);

double XP20(double r, double t, double r0, double D, double v);

double XI20(double r, double t, double r0, double D, double v);

//double XS20(double t, double r0, double D, double v);

double XP30term_nov(double r, double t, double r0, double ka, double D);

double XP30term_v(double r, double t, double r0, double ka, double D, double v);

double XP30(double r, double t, double r0, double ka, double D, double v);

double XI30term_nov(double r, double t, double r0, double ka, double D);

double XI30(double r, double t, double r0, double ka, double D, double v);

double XS30(double t, double r0, double ka, double D, double v);

double XP030(double r, double t, double r0, double ka, double D);

double XI030(double r, double t, double r0, double ka, double D);

double XS030(double t, double r0, double ka, double D);

/* 3D functions below. */

double p_irr(double r, double t, double r0, double kf, double D, double sigma);

double __p_irr(double r, double t, double r0, double D, double sigma, double alpha);

double p_free(double r, double r0, double theta, double t);

double p_survival_irr(double t, double r0, double kf, double D, double sigma);

double __p_reaction_irr(double t, double r0, double kf, double D, double sigma, double alpha, double kD);

double __p_reaction_irr_t_inf(double r0, double kf, double sigma, double kD);

double p_survival_nocollision(double t, double r0, double D, double a);

double dp_survival_nocollision(double t, double r0, double D, double a);

double p_theta_free(double theta, double r, double r0, double t, double D);

double ip_theta_free(double theta, double r, double r0, double t, double D);

/* Functions used in old Brownian Dynamic scheme. */

double GF_EXPORT g_bd_3D(double r0, double sigma, double t, double D);

double GF_EXPORT I_bd_3D(double sigma, double t, double D);

double GF_EXPORT I_bd_r_3D(double r, double sigma, double t, double D);

double GF_EXPORT drawR_gbd_3D(double rnd, double sigma, double t, double D);

double GF_EXPORT g_bd_1D(double r0, double sigma, double t, double D, double v);

double GF_EXPORT I_bd_1D(double sigma, double t, double D, double v);

double GF_EXPORT I_bd_r_1D(double r, double sigma, double t, double D, double v);

double GF_EXPORT drawR_gbd_1D(double rnd, double sigma, double t, double D, double v);

// --------------------------------------------------------------------------------------------------------------------------------
