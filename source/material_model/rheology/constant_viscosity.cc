/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/constant_viscosity.h>
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
      template <int dim>
      ConstantViscosity<dim>::ConstantViscosity ()
      {}


      template <int dim>
      double
      ConstantViscosity<dim>::compute_viscosity (const unsigned int composition) const
      {
        return viscosities[composition];
      }



      template <int dim>
      std::pair<double, double>
      ConstantViscosity<dim>::compute_strain_rate_and_derivative (const double stress,
                                                                  const unsigned int composition) const
      {
        return std::make_pair(0.5*stress/viscosities[composition], 0.5/viscosities[composition]);
      }



      template <int dim>
      void
      ConstantViscosity<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity", "5e24",
                           Patterns::Anything(),
                           "The value of the viscosity $\\eta$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal\\second}.");
      }



      template <int dim>
      void
      ConstantViscosity<dim>::parse_parameters (ParameterHandler &prm,
                                                const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Establish that a background field is required here
        const bool has_background_field = true;

        viscosities = Utilities::parse_map_to_double_array(prm.get("Viscosity"),
                                                           list_of_composition_names,
                                                           has_background_field,
                                                           "Viscosity",
                                                           true,
                                                           expected_n_phases_per_composition);
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
    template class ConstantViscosity<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
