/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/rheology/composite_visco_plastic.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      // The following function is not yet used by ASPECT, so
      // it is commented out to avoid tester issues.
      //
      // It will provide useful output when a material model which
      // uses the CompositeViscoPlastic rheology has been implemented.
      //
      // The composite visco plastic rheology calculates the decomposed strain
      // rates for each of the following deformation mechanisms:
      // diffusion creep, dislocation creep, Peierls creep,
      // Drucker-Prager plasticity and a constant (high) viscosity limiter.
      // The values are provided in this order as a vector of additional
      // outputs. If the user declares one or more mechanisms inactive
      // (by assigning use_mechanism = False) then the corresponding
      // strain rate output will be equal to zero.
      //namespace
      //{
      //  std::vector<std::string> make_strain_rate_additional_outputs_names()
      //  {
      //    std::vector<std::string> names;
      //    names.emplace_back("edot_diffusion");
      //    names.emplace_back("edot_dislocation");
      //    names.emplace_back("edot_peierls");
      //    names.emplace_back("edot_drucker_prager");
      //    names.emplace_back("edot_elasticity");
      //    names.emplace_back("edot_limiter");
      //    return names;
      //  }
      //}


      template <int dim>
      CompositeViscoPlastic<dim>::CompositeViscoPlastic ()
      {}


      template <int dim>
      void
      CompositeViscoPlastic<dim>::fill_outputs(const MaterialModel::MaterialModelInputs<dim> &in,
                                               const std::vector<std::vector<double>> &volume_fractions,
                                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        std::vector<double> average_shear_moduli(in.n_evaluation_points(), 0.);
        const double scale = max_viscosity/(max_viscosity - min_viscosity);
        std::vector<double> partial_strain_rates;

        if (use_elasticity)
          {
            // Compute the effective shear modulus if elasticity is enabled
            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                average_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions[i],
                                                                           elasticity->get_elastic_shear_moduli(),
                                                                           MaterialUtilities::arithmetic);
              }

            // If elasticity is enabled, the effective strain rate must be modified to
            // account for the stored elastic strain from the previous timestep.
            const std::vector<SymmetricTensor<2,dim>> strain_rates_from_stored_stress = elasticity->strain_rates_from_stored_stress(in, average_shear_moduli);

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                // Rescale the effective strain rate to account for the viscosity limiters
                const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i] + scale * strain_rates_from_stored_stress[i];

                // Compute the effective creep viscosity
                const double eta_effcreep = compute_creep_viscosity (in.pressure[i],
                                                                     in.temperature[i],
                                                                     volume_fractions,
                                                                     average_shear_moduli[i],
                                                                     strain_rate,
                                                                     partial_strain_rates,
                                                                     phase_function_values,
                                                                     n_phases_per_composition);
                out.viscosities[i] = 1./(1./max_viscosity + 1./(eta_lim + eta_effcreep));

                // force and reaction will go into elastic_force_outputs and fill_reaction_outputs
                force = - 2. * eta_effcreep / ((eta_lim + eta_effcreep)/max_viscosity + 1.) * strain_rates_from_stored_stress[i];
                tau_creep = 2. * eta_effcreep * strain_rate / (scale + eta_effcreep/max_viscosity);
                reaction = (tau_creep * (1. - elasticity->elastic_damper_viscosity/elasticity->calculate_elastic_viscosity(average_shear_moduli[i]))
                            + 2. * elasticity->elastic_damper_viscosity * strain_rates_from_stored_stress[i]);

                //elasticity->fill_elastic_force_outputs(in, average_shear_moduli, out);
                //elasticity->fill_reaction_outputs(in, average_shear_moduli, out);
              }
            // Fill the elastic additional outputs
            if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
              {
                elastic_out->elastic_shear_moduli = average_shear_moduli;
              }
          }
        else
          {
            // If elasticity is not enabled, the viscosities can be calculated
            // directly from the total strain rate.
            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                // Compute the effective creep viscosity
                const double eta_effcreep = compute_creep_viscosity(in.pressure[i],
                                                                    in.temperature[i],
                                                                    volume_fractions,
                                                                    average_shear_moduli[i],
                                                                    in.strain_rate[i],
                                                                    partial_strain_rates,
                                                                    phase_function_values,
                                                                    n_phases_per_composition);
                out.viscosities[i] = 1./(1./max_viscosity + 1./(eta_lim + eta_effcreep));
              }
          }
      }

      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_creep_viscosity (const double pressure,
                                                           const double temperature,
                                                           const std::vector<double> &volume_fractions,
                                                           const double average_shear_modulus,
                                                           const SymmetricTensor<2,dim> &strain_rate,
                                                           std::vector<double> &partial_strain_rates,
                                                           const std::vector<double> &phase_function_values,
                                                           const std::vector<unsigned int> &n_phases_per_composition) const
      {
        double viscosity = 0.;
        partial_strain_rates.resize(6, 0.);

        // Isostress averaging or the first step toward VRH averaging
        if (viscosity_averaging_scheme == isostress || viscosity_averaging_scheme == voigt_reuss_hill)
          {
            viscosity += compute_isostress_creep_viscosity (pressure,
                                                            temperature,
                                                            volume_fractions,
                                                            average_shear_modulus,
                                                            strain_rate,
                                                            partial_strain_rates,
                                                            phase_function_values,
                                                            n_phases_per_composition);
          }

        // Isostrain averaging or the second step toward VRH averaging
        if (viscosity_averaging_scheme == isostrain || viscosity_averaging_scheme == voigt_reuss_hill)
          {
            double total_volume_fraction = 1.;
            for (unsigned int composition=0; composition < number_of_compositions; ++composition)
              {
                // Only include the contribution to the viscosity
                // from a given composition if the volume fraction exceeds
                // a certain (small) fraction.
                if (volume_fractions[composition] > 1.e-5)
                  {
                    std::vector<double> partial_strain_rates_composition(6, 0.);
                    viscosity += (volume_fractions[composition]
                                  * compute_composition_creep_viscosity (pressure,
                                                                         temperature,
                                                                         composition,
                                                                         average_shear_modulus,
                                                                         strain_rate,
                                                                         partial_strain_rates_composition,
                                                                         phase_function_values,
                                                                         n_phases_per_composition));
                    for (unsigned int j=0; j < 6; ++j)
                      partial_strain_rates[j] += volume_fractions[composition] * partial_strain_rates_composition[j];
                  }
                else
                  {
                    total_volume_fraction -= volume_fractions[composition];
                  }
              }
            viscosity /= total_volume_fraction;
            for (unsigned int j=0; j < 6; ++j)
              partial_strain_rates[j] /= total_volume_fraction;
          }

        // Final step for Voigt-Reuss-Hill averaging is to divide
        // the sum of the viscosities and partial strain rates by two.
        if (viscosity_averaging_scheme == voigt_reuss_hill)
          {
            viscosity /= 2.;
            for (unsigned int j=0; j < 6; ++j)
              partial_strain_rates[j] /= 2.;
          }

        return viscosity;
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_isostress_creep_viscosity (const double pressure,
                                                                     const double temperature,
                                                                     const std::vector<double> &volume_fractions,
                                                                     const double average_shear_modulus,
                                                                     const SymmetricTensor<2,dim> &strain_rate,
                                                                     std::vector<double> &partial_strain_rates,
                                                                     const std::vector<double> &phase_function_values,
                                                                     const std::vector<unsigned int> &n_phases_per_composition) const
      {
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate);

        std::vector<Rheology::DiffusionCreepParameters> diffusion_creep_parameters;
        std::vector<Rheology::DislocationCreepParameters> dislocation_creep_parameters;
        std::vector<Rheology::PeierlsCreepParameters> peierls_creep_parameters;

        // Compute a starting guess at an overall viscosity by taking a
        // harmonic average of the viscosities for each composition.
        // It is assumed that each composition experiences a fraction
        // of the total strain rate corresponding to its volume fraction
        double inveta = 0.;
        for (unsigned int composition=0; composition < number_of_compositions; ++composition)
          {
            const double partial_edot_ii = volume_fractions[composition]*edot_ii;
            double inveta_c = 0.;
            if (use_diffusion_creep)
              {
                diffusion_creep_parameters.push_back(diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition));
                inveta_c += 1./diffusion_creep->compute_viscosity(pressure, temperature, composition, phase_function_values, n_phases_per_composition);
              }

            if (use_dislocation_creep)
              {
                dislocation_creep_parameters.push_back(dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition));
                inveta_c += 1./dislocation_creep->compute_viscosity(partial_edot_ii, pressure, temperature, composition, phase_function_values, n_phases_per_composition);
              }

            if (use_peierls_creep)
              {
                peierls_creep_parameters.push_back(peierls_creep->compute_creep_parameters(composition));
                inveta_c += 1./peierls_creep->compute_approximate_viscosity(partial_edot_ii, pressure, temperature, composition);
              }

            // Crude modification of the creep stress to be no higher than the
            // Drucker-Prager yield stress. Probably fine for a first guess.
            if (use_drucker_prager)
              {
                double creep_stress = 2.*partial_edot_ii/inveta_c;
                const double yield_stress = drucker_prager->compute_yield_stress(drucker_prager_parameters.cohesions[composition],
                                                                                 drucker_prager_parameters.angles_internal_friction[composition],
                                                                                 pressure,
                                                                                 drucker_prager_parameters.max_yield_stress);
                inveta_c = 2.*edot_ii/std::min(creep_stress, yield_stress);
              }

            inveta += inveta_c;
          }
        // First guess at a stress using diffusion, dislocation, and Peierls creep viscosities calculated with the total second strain rate invariant.
        const double eta_guess = std::min(std::max(min_viscosity, 1./inveta), max_viscosity);
        double creep_stress = 2.*eta_guess*edot_ii;

        // In this rheology model, the total strain rate is partitioned between
        // different flow laws. We do not know how the strain is partitioned
        // between these flow laws, nor do we know the creep stress, which is
        // required to calculate the partitioning.

        // The following while loop conducts a Newton iteration to obtain the
        // creep stress, which we need in order to calculate the viscosity.
        double strain_rate_residual = 2*strain_rate_residual_threshold;
        double strain_rate_deriv = 0;
        unsigned int stress_iteration = 0;
        while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {
            double total_volume_fraction = 1.;
            std::pair<double, double> creep_edot_and_deriv = std::make_pair(0., 0.);
            for (unsigned int composition=0; composition < number_of_compositions; ++composition)
              {
                // Only include the contribution to the viscosity
                // from a given composition if the volume fraction exceeds
                // a certain (small) fraction.
                if (volume_fractions[composition] > 1.e-5)
                  {
                    creep_edot_and_deriv = (creep_edot_and_deriv
                                            + volume_fractions[composition]
                                            * compute_creep_strain_rate_and_derivative (creep_stress,
                                                                                        pressure,
                                                                                        temperature,
                                                                                        composition,
                                                                                        diffusion_creep_parameters[composition],
                                                                                        dislocation_creep_parameters[composition],
                                                                                        peierls_creep_parameters[composition],
                                                                                        drucker_prager_parameters,
                                                                                        average_shear_modulus));
                  }
                else
                  {
                    total_volume_fraction -= volume_fractions[composition];
                  }
              }

            creep_edot_and_deriv = std::make_pair(creep_edot_and_deriv.first/total_volume_fraction,
                                                  creep_edot_and_deriv.second/total_volume_fraction);

            const double strain_rate = creep_stress/(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.first;
            strain_rate_deriv = 1./(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.second;

            strain_rate_residual = strain_rate - edot_ii;

            // If the strain rate derivative is zero, we catch it below.
            if (strain_rate_deriv>std::numeric_limits<double>::min())
              creep_stress -= strain_rate_residual/strain_rate_deriv;
            stress_iteration += 1;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(creep_stress)
                                                || !numbers::is_finite(strain_rate_residual)
                                                || !numbers::is_finite(strain_rate_deriv)
                                                || strain_rate_deriv < std::numeric_limits<double>::min()
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " + Utilities::to_string(strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // The creep stress is not the total stress, so we still need to do a little work to obtain the effective viscosity.
        // First, we compute the stress running through the strain rate limiter, and then add that to the creep stress
        // NOTE: The viscosity of the strain rate limiter is equal to (min_visc*max_visc)/(max_visc - min_visc)
        const double lim_stress = 2.*min_viscosity*(edot_ii - creep_stress/(2.*max_viscosity));
        const double total_stress = creep_stress + lim_stress;

        // Compute the strain rate experienced by the different mechanisms
        // These should sum to the total strain rate

        // The components of partial_strain_rates must be provided in the order
        // dictated by make_strain_rate_additional_outputs_names
        double total_volume_fraction = 1.;
        std::fill(partial_strain_rates.begin(), partial_strain_rates.end(), 0);
        for (unsigned int composition=0; composition < number_of_compositions; ++composition)
          {
            // Only include the contribution to the viscosity
            // from a given composition if the volume fraction exceeds
            // a certain (small) fraction.
            if (volume_fractions[composition] > 1.e-5)
              {
                if (use_diffusion_creep)
                  {
                    const std::pair<double, double> diff_edot_and_deriv = diffusion_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, diffusion_creep_parameters[composition]);
                    partial_strain_rates[0] += volume_fractions[composition] * diff_edot_and_deriv.first;
                  }

                if (use_dislocation_creep)
                  {
                    const std::pair<double, double> disl_edot_and_deriv = dislocation_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, dislocation_creep_parameters[composition]);
                    partial_strain_rates[1] += volume_fractions[composition] * disl_edot_and_deriv.first;
                  }

                if (use_peierls_creep)
                  {
                    const std::pair<double, double> prls_edot_and_deriv = peierls_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, peierls_creep_parameters[composition]);
                    partial_strain_rates[2] += volume_fractions[composition] * prls_edot_and_deriv.first;
                  }

                if (use_drucker_prager)
                  {
                    const std::pair<double, double> drpr_edot_and_deriv = drucker_prager->compute_strain_rate_and_derivative(creep_stress, pressure, composition, drucker_prager_parameters);
                    partial_strain_rates[3] += volume_fractions[composition] * drpr_edot_and_deriv.first;
                  }

                if (use_elasticity)
                  {
                    const std::pair<double, double> el_edot_and_deriv = elasticity->compute_strain_rate_and_derivative(creep_stress, average_shear_modulus);
                    partial_strain_rates[4] += volume_fractions[composition] * el_edot_and_deriv.first;
                  }
              }
            else
              {
                total_volume_fraction -= volume_fractions[composition];
              }
          }

        for (unsigned int j=0; j < 4; ++j)
          partial_strain_rates[j] /= total_volume_fraction;

        partial_strain_rates[5] = total_stress/(2.*max_viscosity);

        // Now we return the viscosity using the total stress
        // return total_stress/(2.*edot_ii);

        // Now we return the creep viscosity using the creep stress
        // and creep strain
        return creep_stress/(2.*creep_edot_and_deriv.first);
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_composition_creep_viscosity (const double pressure,
                                                                       const double temperature,
                                                                       const unsigned int composition,
                                                                       const double average_shear_modulus,
                                                                       const SymmetricTensor<2,dim> &strain_rate,
                                                                       std::vector<double> &partial_strain_rates,
                                                                       const std::vector<double> &phase_function_values,
                                                                       const std::vector<unsigned int> &n_phases_per_composition) const
      {
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate);

        Rheology::DiffusionCreepParameters diffusion_creep_parameters;
        Rheology::DislocationCreepParameters dislocation_creep_parameters;
        Rheology::PeierlsCreepParameters peierls_creep_parameters;
        double eta_diff = max_viscosity;
        double eta_disl = max_viscosity;
        double eta_prls = max_viscosity;

        if (use_diffusion_creep)
          {
            diffusion_creep_parameters = diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition);
            eta_diff = diffusion_creep->compute_viscosity(pressure, temperature, composition, phase_function_values, n_phases_per_composition);
          }

        if (use_dislocation_creep)
          {
            dislocation_creep_parameters = dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phases_per_composition);
            eta_disl = dislocation_creep->compute_viscosity(edot_ii, pressure, temperature, composition, phase_function_values, n_phases_per_composition);
          }

        if (use_peierls_creep)
          {
            peierls_creep_parameters = peierls_creep->compute_creep_parameters(composition);
            eta_prls = peierls_creep->compute_approximate_viscosity(edot_ii, pressure, temperature, composition);
          }
        // First guess at a stress using diffusion, dislocation, and Peierls creep viscosities calculated with the total second strain rate invariant.
        const double eta_guess = std::min(std::max(min_viscosity, eta_diff*eta_disl*eta_prls/(eta_diff*eta_disl + eta_diff*eta_prls + eta_disl*eta_prls)), max_viscosity);

        double creep_stress = 2.*eta_guess*edot_ii;

        // Crude modification of the creep stress to be no higher than the
        // Drucker-Prager yield stress. Probably fine for a first guess.
        if (use_drucker_prager)
          {
            const double yield_stress = drucker_prager->compute_yield_stress(drucker_prager_parameters.cohesions[composition],
                                                                             drucker_prager_parameters.angles_internal_friction[composition],
                                                                             pressure,
                                                                             drucker_prager_parameters.max_yield_stress);
            creep_stress = std::min(creep_stress, yield_stress);
          }

        // In this rheology model, the total strain rate is partitioned between
        // different flow laws. We do not know how the strain is partitioned
        // between these flow laws, nor do we know the creep stress, which is
        // required to calculate the partitioning.

        // The following while loop conducts a Newton iteration to obtain the
        // creep stress, which we need in order to calculate the viscosity.
        double strain_rate_residual = 2*strain_rate_residual_threshold;
        double strain_rate_deriv = 0;
        unsigned int stress_iteration = 0;
        while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {

            const std::pair<double, double> creep_edot_and_deriv = compute_creep_strain_rate_and_derivative (creep_stress,
                                                                   pressure,
                                                                   temperature,
                                                                   composition,
                                                                   diffusion_creep_parameters,
                                                                   dislocation_creep_parameters,
                                                                   peierls_creep_parameters,
                                                                   drucker_prager_parameters,
                                                                   average_shear_modulus);

            const double strain_rate = creep_stress/(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.first;
            strain_rate_deriv = 1./(2.*max_viscosity) + (max_viscosity/(max_viscosity - min_viscosity))*creep_edot_and_deriv.second;

            strain_rate_residual = strain_rate - edot_ii;

            // If the strain rate derivative is zero, we catch it below.
            if (strain_rate_deriv>std::numeric_limits<double>::min())
              creep_stress -= strain_rate_residual/strain_rate_deriv;
            stress_iteration += 1;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(creep_stress)
                                                || !numbers::is_finite(strain_rate_residual)
                                                || !numbers::is_finite(strain_rate_deriv)
                                                || strain_rate_deriv < std::numeric_limits<double>::min()
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " + Utilities::to_string(strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // The creep stress is not the total stress, so we still need to do a little work to obtain the effective viscosity.
        // First, we compute the stress running through the strain rate limiter, and then add that to the creep stress
        // NOTE: The viscosity of the strain rate limiter is equal to (min_visc*max_visc)/(max_visc - min_visc)
        const double lim_stress = 2.*min_viscosity*(edot_ii - creep_stress/(2.*max_viscosity));
        const double total_stress = creep_stress + lim_stress;

        // Compute the strain rate experienced by the different mechanisms
        // These should sum to the total strain rate

        // The components of partial_strain_rates must be provided in the order
        // dictated by make_strain_rate_additional_outputs_names
        if (use_diffusion_creep)
          {
            const std::pair<double, double> diff_edot_and_deriv = diffusion_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, diffusion_creep_parameters);
            partial_strain_rates[0] = diff_edot_and_deriv.first;
          }

        if (use_dislocation_creep)
          {
            const std::pair<double, double> disl_edot_and_deriv = dislocation_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, dislocation_creep_parameters);
            partial_strain_rates[1] = disl_edot_and_deriv.first;
          }

        if (use_peierls_creep)
          {
            const std::pair<double, double> prls_edot_and_deriv = peierls_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, peierls_creep_parameters);
            partial_strain_rates[2] = prls_edot_and_deriv.first;
          }

        if (use_drucker_prager)
          {
            const std::pair<double, double> drpr_edot_and_deriv = drucker_prager->compute_strain_rate_and_derivative(creep_stress, pressure, composition, drucker_prager_parameters);
            partial_strain_rates[3] = drpr_edot_and_deriv.first;
          }

        if (use_elasticity)
          {
            const std::pair<double, double> el_edot_and_deriv = elasticity->compute_strain_rate_and_derivative(creep_stress, average_shear_modulus);
            partial_strain_rates[4] += el_edot_and_deriv.first;
          }


        partial_strain_rates[5] = total_stress/(2.*max_viscosity);

        // Now we return the viscosity using the total stress
        // return total_stress/(2.*edot_ii);

        // Now we return the creep viscosity using the creep stress
        // and creep strain
        return creep_stress/(2.*creep_edot_and_deriv.first);
      }



      // Overload the + operator to act on two pairs of doubles.
      std::pair<double,double> operator+(const std::pair<double,double> &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x.first+y.first, x.second+y.second);
      }



      // Overload the + operator to act on a double and a pair of doubles.
      std::pair<double,double> operator+(const double &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x+y.first, x+y.second);
      }



      // Overload the * operator to act on a double and a pair of doubles.
      std::pair<double,double> operator*(const double &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x*y.first, x*y.second);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::compute_creep_strain_rate_and_derivative (const double creep_stress,
                                                                            const double pressure,
                                                                            const double temperature,
                                                                            const unsigned int composition,
                                                                            const DiffusionCreepParameters diffusion_creep_parameters,
                                                                            const DislocationCreepParameters dislocation_creep_parameters,
                                                                            const PeierlsCreepParameters peierls_creep_parameters,
                                                                            const DruckerPragerParameters drucker_prager_parameters,
                                                                            const double average_shear_modulus) const
      {
        std::pair<double, double> creep_edot_and_deriv = std::make_pair(0., 0.);

        if (use_diffusion_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + diffusion_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, diffusion_creep_parameters);

        if (use_dislocation_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + dislocation_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, dislocation_creep_parameters);

        if (use_peierls_creep)
          creep_edot_and_deriv = creep_edot_and_deriv + peierls_creep->compute_strain_rate_and_derivative(creep_stress, pressure, temperature, peierls_creep_parameters);

        if (use_drucker_prager)
          creep_edot_and_deriv = creep_edot_and_deriv + drucker_prager->compute_strain_rate_and_derivative(creep_stress, pressure, composition, drucker_prager_parameters);

        if (use_elasticity)
          creep_edot_and_deriv = creep_edot_and_deriv + elasticity->compute_strain_rate_and_derivative(creep_stress, average_shear_modulus);

        return creep_edot_and_deriv;
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity averaging scheme", "isostrain",
                           Patterns::Selection("isostrain|isostress|Voigt-Reuss-Hill|unspecified"),
                           "When more than one compositional field is present at a point "
                           "with different viscosities, we need to come up with an average "
                           "viscosity at that point.  Select either isostrain "
                           "(the default option), isostress or Voigt-Reuss-Hill.");

        prm.declare_entry ("Include diffusion creep", "true",
                           Patterns::Bool (),
                           "Whether to include diffusion creep in the rheological formulation.");

        prm.declare_entry ("Include dislocation creep", "true",
                           Patterns::Bool (),
                           "Whether to include dislocation creep in the rheological formulation.");

        prm.declare_entry ("Include Peierls creep", "true",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");

        prm.declare_entry ("Include Drucker Prager plasticity", "true",
                           Patterns::Bool (),
                           "Whether to include Drucker-Prager plasticity in the rheological formulation.");

        prm.declare_entry ("Include elasticity", "true",
                           Patterns::Bool (),
                           "Whether to include elasticity in the rheological formulation.");

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Peierls creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        // Drucker Prager parameters
        Rheology::DruckerPrager<dim>::declare_parameters(prm);

        // Elasticity parameters
        Rheology::Elasticity<dim>::declare_parameters(prm);

        // Some of the parameters below are shared with the subordinate
        // rheology models (diffusion, dislocation, ...),
        // and will already have been declared. This is fine, the deal.II
        // parameter handler allows repeated declarations. The default values
        // below will override any defaults given in the subordinate
        // rheology models. The new defaults will apply both to
        // this rheology model and to the subordinate rheology modules.
        prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");

        // Viscosity iteration parameters
        prm.declare_entry ("Strain rate residual tolerance", "1e-22", Patterns::Double(0.),
                           "Tolerance for correct diffusion/dislocation strain rate ratio.");
        prm.declare_entry ("Maximum creep strain rate iterations", "40", Patterns::Integer(0),
                           "Maximum number of iterations to find the correct "
                           "creep strain rate.");

        // Strain rate and stress limiting parameters
        prm.declare_entry ("Minimum viscosity", "1.e17",
                           Patterns::Double(0.),
                           "Minimum effective viscosity. Units: \\si{\\pascal\\second}.");

        prm.declare_entry ("Maximum viscosity", "1.e28",
                           Patterns::Double(0.),
                           "Maximum effective viscosity. Units: \\si{\\pascal\\second}.");
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::parse_parameters (ParameterHandler &prm,
                                                    const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // A background field is required by the subordinate material models
        number_of_compositions = list_of_composition_names.size() + 1;

        if (prm.get ("Viscosity averaging scheme") == "isostrain")
          viscosity_averaging_scheme = isostrain;
        else if (prm.get ("Viscosity averaging scheme") == "isostress")
          viscosity_averaging_scheme = isostress;
        else if (prm.get ("Viscosity averaging scheme") == "Voigt-Reuss-Hill")
          viscosity_averaging_scheme = voigt_reuss_hill;
        else
          {
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));
          }


        min_strain_rate = prm.get_double("Minimum strain rate");

        // Iteration parameters
        strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum creep strain rate iterations");

        // Read min and max viscosity parameters
        min_viscosity = prm.get_double ("Minimum viscosity");
        max_viscosity = prm.get_double ("Maximum viscosity");

        // Rheological parameters

        // Diffusion creep parameters
        use_diffusion_creep = prm.get_bool ("Include diffusion creep");
        if (use_diffusion_creep)
          {
            diffusion_creep = std_cxx14::make_unique<Rheology::DiffusionCreep<dim>>();
            diffusion_creep->initialize_simulator (this->get_simulator());
            diffusion_creep->parse_parameters(prm, expected_n_phases_per_composition);
          }

        // Dislocation creep parameters
        use_dislocation_creep = prm.get_bool ("Include dislocation creep");
        if (use_dislocation_creep)
          {
            dislocation_creep = std_cxx14::make_unique<Rheology::DislocationCreep<dim>>();
            dislocation_creep->initialize_simulator (this->get_simulator());
            dislocation_creep->parse_parameters(prm, expected_n_phases_per_composition);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep");
        if (use_peierls_creep)
          {
            peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm);
          }

        // Drucker Prager parameters
        use_drucker_prager = prm.get_bool ("Include Drucker Prager plasticity");
        if (use_drucker_prager)
          {
            drucker_prager = std_cxx14::make_unique<Rheology::DruckerPrager<dim>>();
            drucker_prager_parameters = drucker_prager->parse_parameters(number_of_compositions, prm);

            AssertThrow(prm.get_bool("Use plastic damper") && prm.get_double("Plastic damper viscosity") > 0.,
                        ExcMessage("If Drucker-Prager plasticity is included in the rheological formulation, you must use a viscous damper with a positive viscosity."));
          }

        // Elasticity parameters
        use_elasticity = prm.get_bool ("Include elasticity");
        if (use_elasticity)
          {
            elasticity = std_cxx14::make_unique<Rheology::Elasticity<dim>>();
            elasticity->initialize_simulator (this->get_simulator());
            elasticity->parse_parameters(prm);
          }

        AssertThrow(use_diffusion_creep == true || use_dislocation_creep == true || use_peierls_creep == true || use_drucker_prager == true || use_elasticity == true,
                    ExcMessage("You need to include at least one deformation mechanism."));

      }

      template <int dim>
      void
      CompositeViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (use_elasticity)
          elasticity->create_elastic_outputs(out);
      }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class CompositeViscoPlastic<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
