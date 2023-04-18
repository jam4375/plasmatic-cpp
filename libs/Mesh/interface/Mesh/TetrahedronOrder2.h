#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class TetrahedronOrder2 : public Element {
  public:
    TetrahedronOrder2(const std::array<Integer, 10> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes);

    virtual Integer NumNodes() const override { return 10; }

    virtual Integer VTKCellType() const override { return 24; }

    virtual Float ShapeFn(Integer index, const Coord &coord) const override;

    Float ShapeFnDerivative(Integer index, Integer dimension, Float xi, Float eta, Float zeta) const;

    virtual Float ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const override;

    virtual Integer GetNodeIndex(Integer index) const override { return _nodeIndices.at(static_cast<size_t>(index)); }

    virtual Float Integrate(const std::function<Float(const Coord &)> integrand) const override;

    static Float ComputeVolume(const Coord &p0, const Coord &p1, const Coord &p2, const Coord &p3);

    std::array<Float, 3> PhysicalToParentCoords(const Coord &coord) const;

    Coord ParentToPhysicalCoords(const std::array<Float, 3> &parent_coords) const;

  private:
    std::array<Integer, 10> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
