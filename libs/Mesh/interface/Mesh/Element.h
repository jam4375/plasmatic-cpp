#pragma once

#include "Coord.h"
#include "Utility/Utility.h"

#include <Eigen/Dense>

#include <array>
#include <vector>

namespace plasmatic {

class Element {
  public:
    virtual ~Element();

    virtual Integer NumNodes() const = 0;

    virtual Integer VTKCellType() const = 0;

    virtual Float ShapeFn(Integer index, const Coord &coord) const = 0;

    virtual Float ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const = 0;

    virtual Integer GetNodeIndex(Integer index) const = 0;

    virtual Float Integrate(const std::function<Float(const Coord &)> integrand) const = 0;

    virtual Eigen::MatrixXd Integrate(const std::function<Eigen::MatrixXd(const Coord &)> integrand, Integer rows,
                                      Integer cols) const = 0;

  private:
};

} // namespace plasmatic
