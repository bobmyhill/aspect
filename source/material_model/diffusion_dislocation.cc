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

#include <aspect/material_model/diffusion_dislocation.h>

using namespace dealii;

namespace aspect
{
  namespace
  {
    std::vector<double>
    get_vector_double (const std::string &parameter, const unsigned int n_fields, ParameterHandler &prm)
    {
      std::vector<double> parameter_list;
      parameter_list = Utilities::string_to_double(Utilities::split_string_list(prm.get (parameter)));
      if (parameter_list.size() == 1)
        parameter_list.resize(n_fields, parameter_list[0]);

      AssertThrow(parameter_list.size() == n_fields,
                  ExcMessage("Length of "+parameter+" list must be either one, or n_compositional_fields+1"));

      return parameter_list;
    }
  }

  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    DiffusionDislocation<dim>::
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
    double
    DiffusionDislocation<dim>::
    average_value ( const std::vector<double> &composition,
                    const std::vector<double> &parameter_values,
                    const enum averaging_scheme &average_type) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }


    template <int dim>
    std::vector<double>
    DiffusionDislocation<dim>::
    calculate_viscosities ( const std::vector<double> &volume_fractions,
                            const double &pressure,
                            const double &temperature,
                            const SymmetricTensor<2,dim> &strain_rate) const
    {
      // Viscosities
      // In the first time step, strain rate is normally zero
      // Sometimes it is also zero at other times.
      // Let's set it to some very small number
      const double e2inv = (strain_rate.norm() == 0.0
                            ?
                            2.0*std::numeric_limits<double>::min()
                            :
                            std::sqrt(abs(second_invariant(strain_rate))));


      // ---- Find effective viscosities for each of the individual phases
      std::vector<double> composition_viscosities(volume_fractions.size()); // viscosities should have same number of entries as compositional fields
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          // \eta = 0.5 * A_i^{-\frac{1}{n_i}} d^\frac{m}{n} e2inv_i^{\frac{1-n_i}{n_i}} \exp\left(\frac{E_i^* + PV_i^*}{n_iRT}\right)
          // where i corresponds to diffusion or dislocation creep

          double e2inv_diffusion = 0;
          double e2inv_dislocation = e2inv;

          // For diffusion creep, viscosity is normally strain rate independent (n=1) but grain size dependent
          // However, occasionally strain rate laws are published with a strain rate dependence,
          // so here we keep the additional term
          double viscosity_diffusion = 0.5 * std::pow(prefactors_diffusion[j], -1/stress_exponents_diffusion[j]) *
                                       std::pow(e2inv_diffusion, ((1-stress_exponents_diffusion[j])/stress_exponents_diffusion[j])) *
                                       std::pow(grain_size, grain_size_exponents_diffusion[j]/stress_exponents_diffusion[j]) *
                                       std::exp((activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
                                                (stress_exponents_diffusion[j]*constants::gas_constant*temperature));

          // For dislocation creep, viscosity is grain size independent (m=0) but strain rate dependent
          double viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation[j], -1/stress_exponents_dislocation[j]) *
                                         std::pow(e2inv_dislocation, ((1-stress_exponents_dislocation[j])/stress_exponents_dislocation[j])) *
                                         std::exp((activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
                                                  (stress_exponents_dislocation[j]*constants::gas_constant*temperature));


          // This loop iterates to find the relative proportions of diffusion and dislocation strain rate
          double viscosity_dislocation_old=0; // old value, start at zero to iterate at least once.
          unsigned int viscosity_iteration = 0;
          while (std::abs((viscosity_dislocation-viscosity_dislocation_old) / viscosity_dislocation)
                 > viscosity_iteration_threshold
                 && viscosity_iteration < viscosity_max_iteration_number)
            {
              viscosity_dislocation_old=viscosity_dislocation;
              e2inv_dislocation = std::max(2.0*std::numeric_limits<double>::min(), viscosity_diffusion / (viscosity_diffusion + viscosity_dislocation) * e2inv);

              if (std::abs(stress_exponents_diffusion[j] - 1) > 1e-5) // the diffusion viscosity only changes if the stress exponent > 1
                {
                  e2inv_diffusion = e2inv - e2inv_dislocation;
                  viscosity_diffusion = 0.5 * std::pow(prefactors_diffusion[j], -1/stress_exponents_diffusion[j]) *
                                        std::pow(e2inv_diffusion, ((1-stress_exponents_diffusion[j])/stress_exponents_diffusion[j])) *
                                        std::pow(grain_size, grain_size_exponents_diffusion[j]/stress_exponents_diffusion[j]) *
                                        std::exp((activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
                                                 (stress_exponents_diffusion[j]*constants::gas_constant*temperature));
                }

              viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation[j], -1/stress_exponents_dislocation[j]) *
                                      std::pow(e2inv_dislocation, ((1-stress_exponents_dislocation[j])/stress_exponents_dislocation[j])) *
                                      std::exp((activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
                                               (stress_exponents_dislocation[j]*constants::gas_constant*temperature));
              viscosity_iteration++;
            }

          // Now find the effective viscosity, with minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(std::pow((1.0/viscosity_diffusion + 1.0/viscosity_dislocation), -1.0), min_visc), max_visc);
        }
      return composition_viscosities;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);

          // Averaging composition-field dependent properties

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

          // calculate effective viscosity
          if (in.strain_rate.size())
            {
              const std::vector<double> composition_viscosities = calculate_viscosities(volume_fractions, pressure, temperature, in.strain_rate[i]);
              const double veff = average_value(composition, composition_viscosities, viscosity_averaging);
              out.viscosities[i] = veff;
            }

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
    DiffusionDislocation<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    DiffusionDislocation<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
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
    DiffusionDislocation<dim>::
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
    DiffusionDislocation<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("DiffusionDislocation");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::Double(0),
                             "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");

          // Viscosity iteration parameters
          prm.declare_entry ("Relative viscosity tolerance", "1e-12", Patterns::Double(0),
                             "Tolerance for correct diffusion/dislocation viscosity ratio.");
          prm.declare_entry ("Viscosity maximum iterations", "40", Patterns::Integer(0),
                             "Maximum number of iterations to find the correct "
                             "diffusion/dislocation viscosity ratio.");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J / (K * kg)$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $1 / K$");

          // Rheological parameters
          prm.declare_entry ("Grain size", "1.e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          // Diffusion creep parameters
          prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");
          prm.declare_entry ("Stress exponents for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_diffusion$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Grain size exponents for diffusion creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of grain size exponents, $m_diffusion$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

          // Dislocation creep parameters
          prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("DiffusionDislocation");
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

          // Iteration parameters
          viscosity_iteration_threshold=prm.get_double ("Relative viscosity tolerance");
          viscosity_max_iteration_number=prm.get_integer ("Viscosity maximum iterations");

          // Equation of state parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

          // ---- Compositional parameters
          grain_size=prm.get_double("Grain size");
          densities=get_vector_double("Densities", n_fields, prm);
          thermal_expansivities=get_vector_double("Thermal expansivities", n_fields, prm);

          // Rheological parameters
          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          // Rheological parameters
          // Diffusion creep parameters (Stress exponents often but not always 1)
          prefactors_diffusion=get_vector_double("Prefactors for diffusion creep", n_fields, prm);
          stress_exponents_diffusion=get_vector_double("Stress exponents for diffusion creep", n_fields, prm);
          grain_size_exponents_diffusion=get_vector_double("Grain size exponents for diffusion creep", n_fields, prm);
          activation_energies_diffusion=get_vector_double("Activation energies for diffusion creep", n_fields, prm);
          activation_volumes_diffusion=get_vector_double("Activation volumes for diffusion creep", n_fields, prm);
          // Dislocation creep parameters (Note the lack of grain size exponents)
          prefactors_dislocation=get_vector_double("Prefactors for dislocation creep", n_fields, prm);
          stress_exponents_dislocation=get_vector_double("Stress exponents for dislocation creep", n_fields, prm);
          activation_energies_dislocation=get_vector_double("Activation energies for dislocation creep", n_fields, prm);
          activation_volumes_dislocation=get_vector_double("Activation volumes for dislocation creep", n_fields, prm);

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    ASPECT_REGISTER_MATERIAL_MODEL(DiffusionDislocation,
                                   "diffusion dislocation",
                                   " An implementation of a viscous rheology including diffusion"
                                   " and dislocation creep."
                                   " Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   " \n\n"
                                   " \\[v_\\text{eff} = \\left(\\frac{1}{v_\\text{eff}^\\text{diff}}+"
                                   " \\frac{1}{v_\\text{eff}^\\text{dis}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_\\text{i} = 0.5 * A^{-\\frac{1}{n_i}} d^\\frac{m_i}{n_i}"
                                   " \\dot{\\varepsilon_i}^{\\frac{1-n_i}{n_i}}"
                                   " \\exp\\left(\\frac{E_i^* + PV_i^*}{n_iRT}\\right)\\]"
                                   " \n\n"
                                   " where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep,"
                                   " $\\dot{\\varepsilon}$ is the square root of the second invariant of the"
                                   " strain rate tensor, $R$ is the gas constant, $T$ is temperature, "
                                   " and $P$ is pressure."
                                   " $A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents"
                                   " $E_i$ are the activation energies and $V_i$ are the activation volumes."
                                   " \n\n"
                                   " The ratio of diffusion to dislocation strain rate is found by iterating"
                                   " to find the dislocation viscosity which satisfies the above equations."
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/DiffusionDislocation'.")
  }
}
