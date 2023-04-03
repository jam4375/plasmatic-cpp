#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Line : public Element {
  public:
    Line(const std::array<Integer, 2> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes);

    virtual Integer NumNodes() const override { return 2; }

    virtual Integer VTKCellType() const override { return 3; }

    virtual Float ShapeFn(Integer index, const Coord &coord) const override;

    virtual Integer GetNodeInex(Integer index) const override { return _nodeIndices.at(static_cast<size_t>(index)); }

  private:
    std::array<Integer, 2> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
