#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Line : public Element {
  public:
    Line(const std::array<Integer, 2> &node_indices, std::shared_ptr<std::vector<Coord>> nodes);

    virtual Integer NumNodes() const override { return 2; }

    virtual Integer VTKCellType() const override { return 3; }

    Float ShapeFn(Integer index, Float xi) const;

    Float ShapeFnDerivative(Integer index, Integer dimension, Float xi) const;

    virtual Float ShapeFn(Integer index, const Coord &coord) const override;

    virtual Float ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const override;

    virtual Integer GetNodeIndex(Integer index) const override { return _nodeIndices.at(static_cast<size_t>(index)); }

    virtual Float Integrate(const std::function<Float(const Coord &)> integrand) const override;

  private:
    std::array<Integer, 2> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
