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


#include <aspect/traction_boundary_conditions/adiabatic_pressure.h>
#include <aspect/global.h>

namespace aspect
{
  namespace TractionBoundaryConditions
  {

    template <int dim>
    Tensor<1,dim>
    AdiabaticPressure<dim>::
    traction (const Point<dim> &p,
              const Tensor<1,dim> &normal_vector) const
    {
      Tensor<1,dim> traction;
      for (unsigned int d=0; d<dim; ++d)
        traction[d] = normal_vector[d]*this->get_adiabatic_conditions().pressure(p);
      return traction;
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace TractionBoundaryConditions
  {
    ASPECT_REGISTER_TRACTION_BOUNDARY_CONDITIONS(AdiabaticPressure,
                                                 "adiabatic pressure",
                                                 "Implementation of a model in which the boundary "
                                                 "traction is given in terms of the adiabatic pressure "
						 "as described in section "
                                                 "``Boundary traction model|AdiabaticPressure''. "
                                                 "\n\n"
                                                 "No parameters should be passed to this model. ")
  }
}
