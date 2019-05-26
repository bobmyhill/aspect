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
  class InterpolatedUniformGridData : public dealii::Function<dim>
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
               const unsigned int component = 0) const override;

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
      double
      value(const Point<dim> &p,
            const unsigned int component = 0) const override;


    private:
      /**
       * The set of interval endpoints in each of the coordinate directions.
       */
      const std::array<std::pair<double, double>, dim> interval_endpoints;

      /**
       * The number of subintervals in each of the coordinate directions.
       */
      const std::array<unsigned int, dim> n_subintervals;

      /**
       * The data that is to be interpolated.
       */
      const Table<dim, double> data_values;

  };
}

namespace aspect
{
  namespace
     {
       // interpolate a data value from a table where ix denotes
       // the (lower) left endpoint of the interval to interpolate
       // in, and p_unit denotes the point in unit coordinates to do so.
       double
       interpolate(const Table<1, double> &data_values,
                   const TableIndices<1> & ix,
                   const Point<1> &        xi)
       {
         return ((1 - xi[0]) * data_values[ix[0]] +
                 xi[0] * data_values[ix[0] + 1]);
       }

       double
       interpolate(const Table<2, double> &data_values,
                   const TableIndices<2> & ix,
                   const Point<2> &        p_unit)
       {
         return (((1 - p_unit[0]) * data_values[ix[0]][ix[1]] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1]]) *
                   (1 - p_unit[1]) +
                 ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1] +
                  p_unit[0] * data_values[ix[0] + 1][ix[1] + 1]) *
                   p_unit[1]);
       }

       double
       interpolate(const Table<3, double> &data_values,
                   const TableIndices<3> & ix,
                   const Point<3> &        p_unit)
       {
         return ((((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2]] +
                   p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2]]) *
                    (1 - p_unit[1]) +
                  ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2]] +
                   p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2]]) *
                    p_unit[1]) *
                   (1 - p_unit[2]) +
                 (((1 - p_unit[0]) * data_values[ix[0]][ix[1]][ix[2] + 1] +
                   p_unit[0] * data_values[ix[0] + 1][ix[1]][ix[2] + 1]) *
                    (1 - p_unit[1]) +
                  ((1 - p_unit[0]) * data_values[ix[0]][ix[1] + 1][ix[2] + 1] +
                   p_unit[0] * data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1]) *
                    p_unit[1]) *
                   p_unit[2]);
       }


       // Interpolate the gradient of a data value from a table where ix
       // denotes the lower left endpoint of the interval to interpolate
       // in, p_unit denotes the point in unit coordinates, and dx
       // denotes the width of the interval in each dimension.
       Tensor<1, 1>
       gradient_interpolate(const Table<1, double> &data_values,
                            const TableIndices<1> & ix,
                            const Point<1> &        p_unit,
                            const Point<1> &        dx)
       {
         (void)p_unit;
         Tensor<1, 1> grad;
         grad[0] = (data_values[ix[0] + 1] - data_values[ix[0]]) / dx[0];
         return grad;
       }


       Tensor<1, 2>
       gradient_interpolate(const Table<2, double> &data_values,
                            const TableIndices<2> & ix,
                            const Point<2> &        p_unit,
                            const Point<2> &        dx)
       {
         Tensor<1, 2> grad;
         double       u00 = data_values[ix[0]][ix[1]],
                u01       = data_values[ix[0] + 1][ix[1]],
                u10       = data_values[ix[0]][ix[1] + 1],
                u11       = data_values[ix[0] + 1][ix[1] + 1];

         grad[0] =
           ((1 - p_unit[1]) * (u01 - u00) + p_unit[1] * (u11 - u10)) / dx[0];
         grad[1] =
           ((1 - p_unit[0]) * (u10 - u00) + p_unit[0] * (u11 - u01)) / dx[1];
         return grad;
       }


       Tensor<1, 3>
       gradient_interpolate(const Table<3, double> &data_values,
                            const TableIndices<3> & ix,
                            const Point<3> &        p_unit,
                            const Point<3> &        dx)
       {
         Tensor<1, 3> grad;
         double       u000 = data_values[ix[0]][ix[1]][ix[2]],
                u001       = data_values[ix[0] + 1][ix[1]][ix[2]],
                u010       = data_values[ix[0]][ix[1] + 1][ix[2]],
                u100       = data_values[ix[0]][ix[1]][ix[2] + 1],
                u011       = data_values[ix[0] + 1][ix[1] + 1][ix[2]],
                u101       = data_values[ix[0] + 1][ix[1]][ix[2] + 1],
                u110       = data_values[ix[0]][ix[1] + 1][ix[2] + 1],
                u111       = data_values[ix[0] + 1][ix[1] + 1][ix[2] + 1];

         grad[0] =
           ((1 - p_unit[2]) *
              ((1 - p_unit[1]) * (u001 - u000) + p_unit[1] * (u011 - u010)) +
            p_unit[2] *
              ((1 - p_unit[1]) * (u101 - u100) + p_unit[1] * (u111 - u110))) /
           dx[0];
         grad[1] =
           ((1 - p_unit[2]) *
              ((1 - p_unit[0]) * (u010 - u000) + p_unit[0] * (u011 - u001)) +
            p_unit[2] *
              ((1 - p_unit[0]) * (u110 - u100) + p_unit[0] * (u111 - u101))) /
           dx[1];
         grad[2] =
           ((1 - p_unit[1]) *
              ((1 - p_unit[0]) * (u100 - u000) + p_unit[0] * (u101 - u001)) +
            p_unit[1] *
              ((1 - p_unit[0]) * (u110 - u010) + p_unit[0] * (u111 - u011))) /
           dx[2];

         return grad;
       }
   } // namespace

  template <int dim>
  inline
  InterpolatedUniformGridData<dim>::InterpolatedUniformGridData(
    const std::array<std::pair<double, double>, dim> &interval_endpoints,
    const std::array<unsigned int, dim>              &n_subintervals,
    const Table<dim, double>                         &data_values)
    :
    interval_endpoints(interval_endpoints),
    n_subintervals(n_subintervals),
    data_values(data_values)
  {
    for (unsigned int d = 0; d < dim; ++d)
      {
        Assert(n_subintervals[d] >= 1,
               ExcMessage("There needs to be at least one subinterval in each "
                          "coordinate direction."));
        Assert(interval_endpoints[d].first < interval_endpoints[d].second,
               ExcMessage("The interval in each coordinate direction needs "
                          "to have positive size"));
        Assert(data_values.size()[d] == n_subintervals[d] + 1,
               ExcMessage("The data table does not have the correct size."));
      }
  }


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

    // find out where this data point lies, relative to the given
    // subdivision points
    TableIndices<dim> ix;
    for (unsigned int d = 0; d < dim; ++d)
      {
        const double delta_x =
          ((this->interval_endpoints[d].second - this->interval_endpoints[d].first) /
              this->n_subintervals[d]);
        if (p[d] <= this->interval_endpoints[d].first)
          ix[d] = 0;
        else if (p[d] >= this->interval_endpoints[d].second - delta_x)
          ix[d] = this->n_subintervals[d] - 1;
        else
          ix[d] = static_cast<unsigned int>(
                    (p[d] - this->interval_endpoints[d].first) / delta_x);
      }

    // now compute the relative point within the interval/rectangle/box
    // defined by the point coordinates found above. truncate below and
    // above to accommodate points that may lie outside the range
    Point<dim> p_unit;
    Point<dim> delta_x;
    for (unsigned int d = 0; d < dim; ++d)
      {
        delta_x[d] =
          ((this->interval_endpoints[d].second - this->interval_endpoints[d].first) /
              this->n_subintervals[d]);
        p_unit[d] = std::max(std::min((p[d] - this->interval_endpoints[d].first -
                                       ix[d] * delta_x[d]) /
                                      delta_x[d],
                                      1.),
                             0.);
      }

    return gradient_interpolate(this->data_values, ix, p_unit, delta_x);
  }

  template <int dim>
  inline
  double
  InterpolatedUniformGridData<dim>::value(const Point<dim> & p,
                                          const unsigned int component) const
  {
    (void)component;
    Assert(
      component == 0,
      ExcMessage(
        "This is a scalar function object, the component can only be zero."));

    // find out where this data point lies, relative to the given
    // subdivision points
    TableIndices<dim> ix;
    for (unsigned int d = 0; d < dim; ++d)
      {
        const double delta_x =
          ((interval_endpoints[d].second - interval_endpoints[d].first) /
           n_subintervals[d]);
        if (p[d] <= interval_endpoints[d].first)
          ix[d] = 0;
        else if (p[d] >= interval_endpoints[d].second - delta_x)
          ix[d] = n_subintervals[d] - 1;
        else
          ix[d] = static_cast<unsigned int>(
            (p[d] - interval_endpoints[d].first) / delta_x);
      }

    // now compute the relative point within the interval/rectangle/box
    // defined by the point coordinates found above. truncate below and
    // above to accommodate points that may lie outside the range
    Point<dim> p_unit;
    for (unsigned int d = 0; d < dim; ++d)
      {
        const double delta_x =
          ((interval_endpoints[d].second - interval_endpoints[d].first) /
           n_subintervals[d]);
        p_unit[d] = std::max(std::min((p[d] - interval_endpoints[d].first -
                                       ix[d] * delta_x) /
                                        delta_x,
                                      1.),
                             0.);
      }

    return interpolate(data_values, ix, p_unit);
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
