/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#ifndef __aspect__model_peierls_h
#define __aspect__model_peierls_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model based on a viscous rheology including diffusion, 
     * dislocation and Peierls law creep. Peierls law creep becomes dominant under
     * low temperature, high strain rate conditions. The expression describing this 
     * creep is dependent on stress, both as a power law function, and as a 
     * factor affecting the activation enthalpy. 
     *
     * The effective viscosity is defined as the harmonic mean of the three 
     * effective viscosity functions describing diffusion, dislocation and Peierls 
     * creep:
     * \\[v_{eff} = \\left(\\frac{1}{v_{eff}^{diff}}+\\frac{1}{v_{eff}^{dis} + \\frac{1}{v_{eff}^{peierls}}\\right)^{-1}\\]
     * where
     * \\[v_{eff}^{diff} = A_{diff}^-1 \\exp\\left(\frac{E_{diff} + PV_{diff}}{RT}\\right)\\]
     * \\[v_{eff}^{dis} =  A_{dis}^{\\frac{-1}{n_dis}} \\dot{\\varepsilon}^{\frac{1-n}{n}} \\exp\\left(\frac{E_{diff} + PV_{diff}}{n_{dis}RT}\\right)\\]
     * \\[v_{eff}^{peierls} = A_{peierls}^{\\frac{-1}{n_peierls}} \\sigma^{\frac{1-n}} \\exp\\left(\frac{E_{diff} + PV_{diff}}{RT} \\left( 1 - \\left(\\frac{\\sigma}{\\sigma_{peierls}} \\right)^p \\right)^q \\right)\\right)\\]
     * \\sigma = v_{eff} \\dot{\\varepsilon}
     *
     * where $\\dot{\\varepsilon}$ is the second invariant of the strain rate tensor,
     * $\\sigma$ is the stress calculated from the viscosities from the previous
     * timestep, $A_i$ are prefactors for each of the creep laws,
     * $E_i$ are the activation energies, $V_i$ are the activation volumes,
     * $n_i$ are stress exponents for dislocation and Peierls creep,
     * $p$ and $q$ are the two other exponents for Peierls creep
     * $\\rho_m$ is the mantle density, $R$ is the gas constant,
     * $T$ is temperature, and $P$ is pressure.
     *
     * Several model parameters (reference densities, thermal expansivities
     * and rheology parameters) can be defined per-compositional field. 
     * If a list of values is given for the density and thermal expansivity, 
     * the weighted sum of the values based on volume fractions of the 
     * compositional fields is used in their place. For the rheological parameters,
     * averaging is saved for the output viscosities. A range of options for 
     * viscosity averaging are possible. If only one value is given for 
     * any of these parameters, all compositions are assigned the same value. 
     * The first value in the list is the value assigned to "background mantle" 
     * (regions where the sum of the compositional fields is < 1.0).
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Peierls : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        typedef typename aspect::MaterialModel::Interface<dim>::MaterialModelInputs MaterialModelInputs;
        typedef typename aspect::MaterialModel::Interface<dim>::MaterialModelOutputs MaterialModelOutputs;

        virtual void evaluate(const MaterialModelInputs &in, MaterialModelOutputs &out) const;
	
        /**
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
	
        double reference_T;
	double min_strain_rate;
	double min_visc;
	double max_visc;
	double veff_coefficient;
	double ref_visc;
	double thermal_diffusivity;
	double heat_capacity;
	
        /**
         * From multicomponent material model: From a list of compositional
         * fields of length N, we come up with an N+1 length list that which
         * also includes the fraction of ``background mantle''. This list
         * should sum to one, and is interpreted as volume fractions.  If the
         * sum of the compositional_fields is greater than one, we assume that
         * there is no background mantle (i.e., that field value is zero).
         * Otherwise, the difference between the sum of the compositional
         * fields and 1.0 is assumed to be the amount of background mantle.
         */
        std::vector<double> compute_volume_fractions(
	   const std::vector<double> &compositional_fields) const;
        std::vector<double> densities;
        std::vector<double> thermal_expansivities;

        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         * Select between harmonic, arithmetic, geometric, and
         * maximum_composition.  The max composition scheme simply uses the
         * viscosity of whichever field has the highes volume fraction.
         */
        enum averaging_scheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        } viscosity_averaging;


        virtual double average_value (const std::vector<double> &composition,
				      std::vector<double> &parameter_values,
				      const enum averaging_scheme &average_type) const;

	virtual std::vector<double> calculate_viscosities ( const std::vector<double> &volume_fractions,
							    const double &pressure,
							    const double &temperature,
							    const SymmetricTensor<2,dim> &strain_rate,
							    const double &old_viscosity ) const;

	
	std::vector<double> prefactors_diffusion;
	std::vector<double> activation_energies_diffusion;
	std::vector<double> activation_volumes_diffusion;

	std::vector<double> prefactors_dislocation;
	std::vector<double> activation_energies_dislocation;
	std::vector<double> activation_volumes_dislocation;
	std::vector<double> stress_exponents_dislocation;

	std::vector<double> prefactors_Peierls;
	std::vector<double> activation_energies_Peierls;
	std::vector<double> activation_volumes_Peierls;
	std::vector<double> stress_exponents_Peierls;
	std::vector<double> reference_stresses_Peierls;
	std::vector<double> p_exponents_Peierls;
	std::vector<double> q_exponents_Peierls;
    };

  }
}

#endif
