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


#ifndef __aspect__traction_boundary_conditions_adiabatic_pressure_h
#define __aspect__traction_boundary_conditions_adiabatic_pressure_h

#include <aspect/simulator_access.h>
#include <aspect/traction_boundary_conditions/interface.h>

namespace aspect
{
  namespace TractionBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements adiabatic pressure traction boundary 
     * conditions. This is equivalent to an open boundary condition 
     * in domains where hydrostatic pressure is relevant.
     *
     * @ingroup TractionBoundaryConditionsModels
     */
    template <int dim>
    class AdiabaticPressure : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         *
         * For the current class, this function returns a tensor equal
	 * to the adiabatic pressure.
         */
        virtual
        Tensor<1,dim>
        traction (const Point<dim> &position,
                  const Tensor<1,dim> &normal_vector) const;
    };
  }
}


#endif
