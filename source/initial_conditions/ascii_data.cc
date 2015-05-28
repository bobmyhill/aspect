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


#include <aspect/global.h>
#include <aspect/initial_conditions/ascii_data.h>

namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(1);
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      return Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-conditions/ascii-data/test/",
                                                          "box_2d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(AsciiData,
                                       "ascii data",
                                       "Implementation of a model in which the initial temperature "
                                       "is derived from a file containing data in ascii format. "
				       "ASPECT ignores the first lines if they begin with '#', "
				       "except for a required line starting '# POINTS:', "
				       " which contains the number of grid points in each dimension, "
				       " e.g.. '# POINTS: 3 3 3'. "
				       "ASPECT then reads one datum per line, in three/four column format: "
				       "'x', 'y', '(z)', 'Temperature [K]' for Cartesian data, and "
				       "'r', 'phi' [longitude], '(theta)' [pi/2-latitude], 'Temperature [K]' "
				       "for spherical data. Angular data ('phi' and 'theta') should be "
				       "supplied in radians. The 'x' variable is in the innermost loop, "
				       "followed by 'y' (and 'z'). In other words, the first two data lines "
				       "will be x[0] y[0] (z[0]) T[0,0,0] and x[1] y[0] (z[0]) T[1,0,0].")
  }
}
