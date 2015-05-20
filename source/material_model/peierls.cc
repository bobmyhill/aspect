/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/peierls.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    Peierls<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions( compositional_fields.size()+1);

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    void
    Peierls<dim>::
    evaluate(const MaterialModelInputs &in,
             MaterialModelOutputs &out) const
    {
      const double R = 8.31446; // J mol-1 K-1
      const double g = 9.80665; // m s-2
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
	  const double pressure= in.pressure[i] // NEED TO CHECK THIS
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i];
	  const viscosity = in.viscosities[i]; // Try using viscosity from the previous time step to compute stresses for this time step
	  
	  // densities
          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];


	  // prefactors_diffusion
          // activation_energies_diffusion
          // activation_volumes_diffusion
          // prefactors_dislocation
          // activation_energies_dislocation
          // activation_volumes_dislocation
          // stress_exponents_dislocation
          // prefactors_Peierls
          // activation_energies_Peierls
          // activation_volumes_Peierls
          // stress_exponents_Peierls
	  // reference_stresses_Peierls
	  // p_exponents_Peierls
	  // q_exponents_Peierls






	  
          const double e2inv = std::sqrt(std::pow(second_invariant(strain_rate),2) + std::pow(min_strain_rate,2)); // include minimum strain rate
	  const double sigma2inv = viscosity/e2inv; // second invariant of the stress
	  const double non_dimensional_peierls_stress = std::pow(sigma2inv/reference_peierls_stress, peierls_exponent_p); 
	  const double peierls_factor = std::pow(1.0-non_dimensional_peierls_stress, peierls_exponent_q);
	  
	  
	  // Find the individual components of the viscosities
	  const double viscosity_diffusion = std::min(1e22,(1e0/prefactor_diffusion)*
						      std::exp((activation_energy_diffusion+activation_volume_diffusion*pressure)/(R*temperature)));
	  
	  const double viscosity_dislocation = std::min(1e22,std::pow(prefactor_dislocation,-1e0/stress_exponent_dislocation)*
							std::pow(e2inv,(1e0-stress_exponent_dislocation)/
								 stress_exponent_dislocation)*
							std::exp((activation_energy_dislocation+
								  activation_volume_dislocation*pressure)/(stress_exponent_dislocation*R*temperature)));
	  
	  const double viscosity_peierls = std::min(1e22,(1e0/prefactor_peierls)*
						    std::pow(sigma2inv, 1e0-stress_exponent_peierls)*
						    std::exp(peierls_factor*(activation_energy_peierls+activation_volume_peierls*pressure)/(R*temperature)));

											  
          // Effective viscosity = harmonic mean of diffusion, dislocation and Peierls creep. Range is limited to 1e17-1e28 for stability.
          const double veff = std::min(std::max(std::pow((1.0/viscosity_diffusion + 1.0/viscosity_dislocation + 1.0/viscosity_peierls), -1.0), min_visc), max_visc);

          out.viscosities[i] = veff;
          out.densities[i] = density;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = heat_capacity;
          // Thermal conductivity at the given positions.
          out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity * density;
          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    Peierls<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    Peierls<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    Peierls<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the viscosity() function
      // to see the dependencies
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none))
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else if (((dependence & NonlinearDependence::pressure) != NonlinearDependence::none))
        return true;
      else if (((dependence & NonlinearDependence::strain_rate) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Peierls<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else if (((dependence & NonlinearDependence::pressure) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Peierls<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Peierls<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Peierls<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Peierls<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    Peierls<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Peierls");
        {
	  // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0), "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0), "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0), "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0), "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::Double(0), "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0), "Reference viscosity for nondimensionalization.");

	 
	  // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J / (K * kg)$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / K$");

	  // Rheological parameters
	  // Diffusion creep parameters
	  prm.declare_entry ("Prefactors for diffusion creep", "1.92e-11",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / s$");
          prm.declare_entry ("Activation energies for diffusion creep", "335e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
	  prm.declare_entry ("Activation volumes for diffusion creep", "6.4e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

	  // Dislocation creep parameters
	  prm.declare_entry ("Prefactors for dislocation creep", "2.42e-10",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $Pa^{-n_dislocation} s^{-1}$");
	  prm.declare_entry ("Activation energies for dislocation creep", "540e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
	  prm.declare_entry ("Activation volumes for dislocation creep", "6.4e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponents for dislocation creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");

	  // Peierls creep parameters
	  prm.declare_entry ("Prefactors for Peierls creep", "1.00e2",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $Pa^{-n_peierls}s^{-1}");
	  prm.declare_entry ("Activation energies for Peierls creep", "7.66e5",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
	  prm.declare_entry ("Activation volumes for Peierls creep", "0.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponents for Peierls creep", "2",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_peierls$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Reference stresses for Peierls creep", "4.2e9",
			     Patterns::List(Patterns::Double(0)),
			     "($\\sigma_{peierls}$). Units: Pa");
          prm.declare_entry ("p exponents for Peierls creep", "1",
			     Patterns::List(Patterns::Double(0)),
			     "(peierls_exponent_p). Units: None");
	  prm.declare_entry ("q exponents for Peierls creep", "2",
			     Patterns::List(Patterns::Double(0)),
			     "(peierls_exponent_q). Units: None");
	  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Peierls<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Peierls");
        {
	  // Initialise empty vector for compositional field variables
          std::vector<double> x_values;

	  // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          min_strain_rate = prm.get_double("Minimum strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          veff_coefficient = prm.get_double ("Effective viscosity coefficient");
          ref_visc = prm.get_double ("Reference viscosity");

	  
	  // Equation of state parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

	  // ---- Densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of density list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            densities.assign( n_fields , x_values[0]);
          else
            densities = x_values;

          // ---- Thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of thermal expansivity list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            thermal_expansivities.assign( n_fields , x_values[0]);
          else
            thermal_expansivities = x_values;


	  // Rheological parameters
	  // Diffusion creep parameters
	  
          // ---- diffusion creep prefactors
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Prefactors for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of prefactors for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            prefactors_diffusion.assign( n_fields , x_values[0] );
          else
            prefactors_diffusion = x_values;

	  // ---- diffusion creep activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies_diffusion.assign( n_fields , x_values[0] );
          else
            activation_energies_diffusion = x_values;

	  // ---- diffusion creep activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volume for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_volumes_diffusion.assign( n_fields , x_values[0] );
          else
            activation_volumes_diffusion = x_values;


	  // Dislocation creep parameters
          // ---- dislocation creep prefactors
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Prefactors for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of prefactors for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            prefactors_dislocation.assign( n_fields , x_values[0] );
          else
            prefactors_dislocation = x_values;

	  // ---- dislocation creep activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies_dislocation.assign( n_fields , x_values[0] );
          else
            activation_energies_dislocation = x_values;

	  // ---- dislocation creep activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volume for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_volumes_dislocation.assign( n_fields , x_values[0] );
          else
            activation_volumes_dislocation = x_values;

	  // ---- dislocation creep stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of stress exponents for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            stress_exponents_dislocation.assign( n_fields , x_values[0] );
          else
            stress_exponents_dislocation = x_values;


	  // Peierls creep parameters
          // ---- Peierls creep prefactors
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Prefactors for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of prefactors for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            prefactors_Peierls.assign( n_fields , x_values[0] );
          else
            prefactors_Peierls = x_values;

	  // ---- Peierls creep activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies_Peierls.assign( n_fields , x_values[0] );
          else
            activation_energies_Peierls = x_values;

	  // ---- Peierls creep activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volume for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_volumes_Peierls.assign( n_fields , x_values[0] );
          else
            activation_volumes_Peierls = x_values;

	  // ---- Peierls creep stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of stress exponents for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            stress_exponents_Peierls.assign( n_fields , x_values[0] );
          else
            stress_exponents_Peierls = x_values;

	  // ---- Peierls reference stresses
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Reference stresses for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of reference stress for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            reference_stresses_Peierls.assign( n_fields , x_values[0] );
          else
            reference_stresses_Peierls = x_values;

	  // ---- Peierls exponents p
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("p exponents for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of p exponents for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            p_exponents_Peierls.assign( n_fields , x_values[0] );
          else
            p_exponents_Peierls = x_values;

	  // ---- Peierls exponents q
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("q exponents for Peierls creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of q exponents for Peierls list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            q_exponents_Peierls.assign( n_fields , x_values[0] );
          else
            q_exponents_Peierls = x_values;
	  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    
    ASPECT_REGISTER_MATERIAL_MODEL(Peierls,
                                   "Peierls",
                                   "An implementation of a viscous rheology including Peierls creep, which incorporates a"
                                   " stress-dependence of the activation energy, and is particularly relevant at low temperatures."
                                   " Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   "\n\n"
                                   " \\[v_{eff} = \\left(\\frac{1}{v_{eff}^{diff}}+\\frac{1}{v_{eff}^{dis} + \\frac{1}{v_{eff}^{peierls}}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_{eff}^{diff} = A_{diff}^-1 \\exp\\left(\frac{E_{diff} + PV_{diff}}{RT}\\right)\\]"
				   " \\[v_{eff}^{dis} =  A_{dis}^{\\frac{-1}{n_dis}} \\dot{\\varepsilon}^{\frac{1-n}{n}} \\exp\\left(\frac{E_{diff} + PV_{diff}}{n_{dis}RT}\\right)\\]"
				   " \\[v_{eff}^{peierls} = A_{peierls}^{\\frac{-1}{n_peierls}} \\dot{\\sigma}^{\frac{1-n}} \\exp\\left(\frac{E_{diff} + PV_{diff}}{RT} \\left( 1 - \\left(\\frac{\\sigma}{\\sigma_{peierls}} \\right)^p \\right)^q \\right)\\right)\\]"
                                   "\n\n"
                                   " Where $\\dot{\\varepsilon}$ is related to the second invariant of the strain rate tensor,"
				   " $\\dot{\\varepsilon}_{ref}$ is a reference strain rate, $n_v$ and $n_p$ are stress exponents,"
				   " $E_a$ is the activation energy, $V_a$ is the activation volume,"
				   " $\\rho_m$ is the mantle density, $R$ is the gas constant,"
				   " $T$ is temperature, and $P$ is pressure."
                                   " \n\n"
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/Peierls'.")
  }
}
