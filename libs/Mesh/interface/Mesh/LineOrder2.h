#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class LineOrder2 : public Element {
  public:
    LineOrder2(const std::array<Integer, 3> &node_indices, std::shared_ptr<std::vector<Coord>> nodes);

    virtual Integer NumNodes() const override { return 3; }

    virtual Integer VTKCellType() const override { return 21; }

    Float ShapeFnDerivative(Integer index, Integer dimension, Float xi) const;

    virtual Float ShapeFn(Integer index, const Coord &coord) const override;

    virtual Float ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const override;

    virtual Integer GetNodeIndex(Integer index) const override { return _nodeIndices.at(static_cast<size_t>(index)); }

    virtual Float Integrate(const std::function<Float(const Coord &)> integrand) const override;

    static Float ComputeLength(const Coord &p0, const Coord &p1);

    Float PhysicalToParentCoords(const Coord &coord) const;

    Coord ParentToPhysicalCoords(const Float &xi) const;

  private:
    std::array<Integer, 3> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
