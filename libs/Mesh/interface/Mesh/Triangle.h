#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Triangle : public Element {
  public:
    Triangle(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes);

    virtual Integer NumNodes() const override { return 3; }

    virtual Integer VTKCellType() const override { return 5; }

    virtual Float ShapeFn(Integer index, const Coord &coord) const override;

    Float ShapeFnDerivative(Integer index, Integer dimension, Float lambda1, Float lambda2) const;

    virtual Float ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const override;

    virtual Integer GetNodeIndex(Integer index) const override { return _nodeIndices.at(static_cast<size_t>(index)); }

    virtual Float Integrate(const std::function<Float(const Coord &)> integrand) const override;

    static Float ComputeArea(const Coord &p1, const Coord &p2, const Coord &p3);

    std::array<Float, 3> PhysicalToParentCoords(const Coord &coord) const;

  private:
    std::array<Integer, 3> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
