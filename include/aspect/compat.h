/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>

// C++11 related includes.
#include <array>
#include <functional>
#include <memory>

#if !DEAL_II_VERSION_GTE(9,2,0)

#include <deal.II/base/table.h>
#include <deal.II/base/function_lib.h>
namespace aspect
{
  using namespace dealii;
  /**
   * A scalar function that computes its values by (bi-, tri-)linear
   * interpolation from a set of point data that are arranged on a uniformly
   * spaced tensor product mesh. This function is derived from
   * deal.II (dealii/include/deal.II/base/function_lib.h)
   */
  template <int dim>
  class InterpolatedUniformGridData : public dealii::Functions::InterpolatedUniformGridData<dim>
  {
    public:
      /**
       * Constructor
       * @param interval_endpoints The left and right end points of the
       * (uniformly subdivided) intervals in each of the coordinate directions.
       * @param n_subintervals The number of subintervals in each coordinate
       * direction. A value of one for a coordinate means that the interval is
       * considered as one subinterval consisting of the entire range. A value
       * of two means that there are two subintervals each with one half of the
       * range, etc.
       * @param data_values A dim-dimensional table of data at each of the mesh
       * points defined by the coordinate arrays above. Note that the Table
       * class has a number of conversion constructors that allow converting
       * other data types into a table where you specify this argument.
       */
      InterpolatedUniformGridData(
        const std::array<std::pair<double, double>, dim> &interval_endpoints,
        const std::array<unsigned int, dim>              &n_subintervals,
        const Table<dim, double>                         &data_values);

      /**
       * Compute the value of the function set by bilinear interpolation of the
       * given data set.
       *
       * @param p The point at which the function is to be evaluated.
       * @param component The vector component. Since this function is scalar,
       * only zero is a valid argument here.
       * @return The interpolated value at this point. If the point lies outside
       * the set of coordinates, the function is extended by a constant.
       */
      virtual
      Tensor<1, dim>
      gradient(const Point<dim> &p,
               const unsigned int component) const;

  };
}

namespace aspect
{
  template <int dim>
  inline
  InterpolatedUniformGridData<dim>::InterpolatedUniformGridData(
    const std::array<std::pair<double, double>, dim> &interval_endpoints,
    const std::array<unsigned int, dim>              &n_subintervals,
    const Table<dim, double>                         &data_values)
    :
    dealii::Functions::InterpolatedUniformGridData<dim>(&interval_endpoints, &n_subintervals, &data_values)
  {}


  /**
   * This function is derived from
   * deal.II (dealii/source/base/function_lib.cc)
  */
  template <int dim>
  inline
  Tensor<1, dim>
  InterpolatedUniformGridData<dim>::gradient(
    const Point<dim> &p,
    const unsigned int component) const
  {
    (void)component;
    Assert(
      component == 0,
      ExcMessage(
        "This is a scalar function object, the component can only be zero."));

    // find out where this data point lies
    const TableIndices<dim> ix = table_index_of_point(p);

    Point<dim> dx;
    for (unsigned int d = 0; d < dim; ++d)
      dx[d] = (this->interval_endpoints[d].second - this->interval_endpoints[d].first)/this->n_subintervals[d];

    Point<dim> p_unit;
    for (unsigned int d = 0; d < dim; ++d)
      p_unit[d] =
        std::max(std::min((p[d] - this->coordinate_values[d][ix[d]]) / dx[d], 1.),
                 0.0);

    return gradient_interpolate(this->data_values, ix, p_unit, dx);
  }

}
#endif


// We would like to use a function from SolverControl that was introduced after
// deal.II 8.5. For older versions use this derived class instead that implements
// the function and the typedef guarantees it is found before deal.II's class.
#if !DEAL_II_VERSION_GTE(9,0,0)

#include <deal.II/lac/solver_control.h>

namespace aspect
{
  using namespace dealii;

  class SolverControl : public dealii::SolverControl
  {
    public:
      SolverControl(const unsigned int n           = 100,
                    const double       tol         = 1.e-10,
                    const bool         log_history = false,
                    const bool         log_result  = true)
        :
        dealii::SolverControl (n, tol, log_history, log_result)
      {}

      dealii::SolverControl::State
      check (const unsigned int step,
             const double check_value)
      {
        dealii::SolverControl::State return_value = dealii::SolverControl::check(step, check_value);

        if (step == 0)
          history_data.resize(history_data.size()+1);

        return return_value;
      }


      const std::vector<double> &get_history_data() const
      {
        Assert (history_data_enabled, ExcHistoryDataRequired());
        Assert (history_data.size() > 0,
                ExcMessage("The SolverControl object was asked for the solver history "
                           "data, but there is no data. Possibly you requested the data before the "
                           "solver was run."));

        return history_data;
      }
  };
}

#include <deal.II/grid/tria_boundary_lib.h>

/**
 * The ConeBoundary implements the normal_vector function.
 */
namespace aspect
{
  using namespace dealii;

  template <int dim>
  class ConeBoundary : public dealii::ConeBoundary<dim>
  {
    public:

      /**
       * Constructor. Here the boundary object is constructed. The points
       * <tt>x_0</tt> and <tt>x_1</tt> describe the starting and ending points of
       * the axis of the (truncated) cone. <tt>radius_0</tt> denotes the radius
       * corresponding to <tt>x_0</tt> and <tt>radius_1</tt> the one corresponding
       * to <tt>x_1</tt>.
       */
      ConeBoundary (const double radius_0,
                    const double radius_1,
                    const Point<dim> x_0,
                    const Point<dim> x_1);

      virtual
      Tensor<1,dim>
      normal_vector (const typename Triangulation<dim>::face_iterator &face,
                     const Point<dim> &p) const;

  };

}
namespace aspect
{
  template <int dim>
  inline
  ConeBoundary<dim>::ConeBoundary (const double radius_0,
                                   const double radius_1,
                                   const Point<dim> x_0,
                                   const Point<dim> x_1)
    :
    dealii::ConeBoundary<dim>(radius_0, radius_1, x_0, x_1)
  {}

  template <int dim>
  inline
  Tensor<1,dim>
  ConeBoundary<dim>::
  normal_vector (const typename Triangulation<dim>::face_iterator &,
                 const Point<dim> &p) const
  {
    // TODO only for cone opening along z-axis
    AssertThrow (dim == 3, ExcInternalError());
    AssertThrow (this->radius_0 == 0., ExcInternalError());
    AssertThrow (this->x_0[0] == 0., ExcInternalError());
    const double c_squared = (this->radius_1 / this->x_1[dim-1])*(this->radius_1 / this->x_1[dim-1]);
    Tensor<1,dim> normal = p;
    normal[0] *= -2.0/c_squared;
    normal[1] *= -2.0/c_squared;
    normal[dim-1] *= 2.0;

    return normal/normal.norm();
  }
}
#endif

#if !DEAL_II_VERSION_GTE(9,0,0)
namespace dealii
{
  namespace std_cxx14
  {
#ifdef DEAL_II_WITH_CXX14
    using std::make_unique;
#else
    /**
     * An implementation of std::make_unique(). Since we require C++11
     * even if the compiler does not support C++14, we can implement
     * everything that they forgot to put into C++11.
     */
    template <typename T, class... Args>
    inline
    std::unique_ptr<T>
    make_unique(Args &&... args)
    {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
#endif
  }
}
#else
#include <deal.II/base/std_cxx14/memory.h>
#endif


#endif
