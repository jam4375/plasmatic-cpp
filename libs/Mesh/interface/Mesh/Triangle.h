#pragma once

#include "Element.h"
#include "Utility/Utility.h"

#include <array>
#include <vector>

namespace plasmatic {

class Triangle : public Element {
  public:
    Triangle(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes);

    virtual Float ShapeFn(Integer index, const Coord &coord) override;

  private:
    std::array<Integer, 3> _nodeIndices;
    std::shared_ptr<std::vector<Coord>> _nodes;
};

} // namespace plasmatic
