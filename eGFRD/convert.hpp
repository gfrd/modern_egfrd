#ifndef CONVERT_HPP
#define CONVERT_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#define N_A          (6.022140857E23)                 // Avogadros number

// --------------------------------------------------------------------------------------------------------------------------------

class convert
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------

   // Convert a reaction rate from units 'per molar per second' to units 'meter^3 per second'
   static double per_M_per_sec_to_m3_per_sec(double rate) { return rate / (1E3 * N_A); }

   // Convert a reaction rate from units 'per millimolar per second' to units 'meter^3 per second'
   static double per_mM_per_sec_to_m3_per_sec(double rate) { return per_M_per_sec_to_m3_per_sec(rate * 1E3); }

   // Convert a reaction rate from units 'per micromolar per second' to units 'meter^3 per second'
   static double per_uM_per_sec_to_m3_per_sec(double rate) { return per_M_per_sec_to_m3_per_sec(rate * 1E6); }

   // Convert a reaction rate from units 'per nanomolar per second' to units 'meter^3 per second'
   static double per_nM_per_sec_to_m3_per_sec(double rate) { return per_M_per_sec_to_m3_per_sec(rate * 1E9); }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Convert a concentration from units 'molar' to units 'per meter^3'
   static double M_to_per_m3(double molar) { return molar * 1E3 * N_A; }

   // Convert a concentration from units 'millimolar' to units 'per meter^3'
   static double mM_to_per_m3(double molar) { return M_to_per_m3(molar * 1E-3); }

   // Convert a concentration from units 'micromolar' to units 'per meter^3'
   static double uM_to_per_m3(double molar) { return M_to_per_m3(molar * 1E-6); }

   // Convert a concentration from units 'nanomolar' to units 'per meter^3'
   static double nM_to_per_m3(double molar) { return M_to_per_m3(molar * 1E-9); }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Calculate the number of particles in a volume (in liter or dm^3) with a concentration (in molar or mol/dm3).
   static double particles_in_volume_liter(double concentration, double volume) { return concentration * volume * N_A; }

   // Calculate the number of particles in a volume (in m^3) with a concentration (in molar or mol/dm3).
   static double particles_in_volume_m3(double concentration, double volume) { return concentration * volume * 1E3 * N_A; }

   // Calculate the number of particles in a volume (in um^3) with a concentration (in molar or mol/dm3).
   static double particles_in_volume_um3(double concentration, double volume) { return concentration * volume * 1E-15 * N_A; }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Calculate the 'pseudo-'reaction rate (kD) caused by diffusion
   // kD is equal to 1 divided by the time it takes for two particles to meet each other by diffusion. It is needed when converting from
   //    an intrinsic reaction rate to an overall reaction rates or vice versa.
   // Example: A + B-> C
   // Arguments
   //   - Dtot : the diffusion constant of particle A plus the diffusion constant of particle B.  [Units: meter^2 / second]
   //   - sigma : the radius of particle A plus the radius of particle B. [Units: meter]
   // Note: This function is only available for reaction rules in 3D.

   static double k_D(double Dtot, double sigma) { return 4.0 * M_PI * Dtot * sigma; }

   // --------------------------------------------------------------------------------------------------------------------------------

   /* Convert an overall reaction rate (kon) for a binding/annihilation reaction rule to an intrinsic reaction rate (ka)

      Example:
      - A + B-> C              binding reaction rule
      - A + B -> 0             annihilation reaction rule

      Arguments:
       - kon                   the overall reaction rate for the reaction rule.  [Units: meters^3 / second]
       - kD                    the 'pseudo-'reaction rate caused by the diffusion of particles A and B. See the function k_D(). [Units: meters^3 / second]

      Note: This function is only available for reaction rules in 3D. No analytical expression for kD in 1D or 2D is currently known.
      */
   static double k_a(double kon, double kD) { if (kon > kD) throw std::runtime_error("kon > kD.");       return 1.0 / (1.0 / kon - 1.0 / kD); }



   /* Convert an overall reaction rate (koff) for an unbinding reaction rule to an intrinsic reaction rate(kd).

      This one is a bit tricky. We consider reaction rules with only 1 reactant.In case there is only 1 product also, no conversion in
      necessary.But when the reaction rule has 2 products, we need to take the reverse reaction rule into account and do the proper conversion.

      Example:
       - C -> A + B            unbinding reaction rule
       - A + B -> C            reverse reaction rule

      Arguments:
       - koff                  the overall reaction rate for the unbinding reaction rule.  [Units: meters^3 / second]
        - kon                  the overall reaction rate for the reverse reaction rule. [Units: meters^3 / second]
        - kD                   the 'pseudo-'reaction rate caused by the diffusion of particles A and B.See the function k_D(). [Units: meters^3 / second]

      Note: This function is only available for reaction rules in 3D. No analytical expression for kD in 1D or 2D is currently known.
      */
   static double k_d(double koff, double kon, double kD) { return k_d_using_ka(koff, k_a(kon, kD), kD); }



   /*Convert an overall reaction rate (koff) for an unbinding reaction
   rule to an intrinsic reaction rate(kd).

   Similar to the function k_d(), but expects an intrinsic rate(ka)
   instead of an overall rate(kon) for the reversed reaction rule as
   the second argument.

   This function is only available for reaction rules in 3D.No
   analytical expression for kD in 1D or 2D is currently known.

   */
   static double k_d_using_ka(double koff, double ka, double kD) { return koff * (1 + ka / kD); }

   /*Convert an intrinsic reaction rate (ka) for a binding/annihilation
   reaction rule to an overall reaction rate(kon).

   The inverse of the function k_a().

   Rarely needed.

   This function is only available for reaction rules in 3D.No
   analytical expression for kD in 1D or 2D is currently known.

   */
   static double k_on(double ka, double kD) { return 1.0 / ((1.0 / kD) + (1.0 / ka)); }

   /*Convert an intrinsic reaction rate (kd) for an unbinding reaction
   rule to an overall reaction rate(koff).

   The inverse of the function k_d().

   Rarely needed.

   This function is only available for reaction rules in 3D.No
   analytical expression for kD in 1D or 2D is currently known.

   */
   static double k_off(double kd, double kon, double kD) { return k_off_using_ka(kd, k_a(kon, kD), kD); }

   /*Convert an intrinsic reaction rate (kd) for an unbinding reaction
   rule to an overall reaction rate(koff).

   Similar to the function k_off(), but expects an intrinsic rate
   (ka) instead of an overall rate(kon) as the second argument.

   Rarely needed.

   This function is only available for reaction rules in 3D.No
   analytical expression for kD in 1D or 2D is currently known.

   */

   static double k_off_using_ka(double kd, double ka, double kD) { return 1.0 / (float(ka) / (kd * kD) + (1.0 / kd)); }

};


// --------------------------------------------------------------------------------------------------------------------------------

#endif /* CONVERT_HPP */
