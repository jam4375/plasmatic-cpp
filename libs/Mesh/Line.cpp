#include "interface/Mesh/Line.h"

namespace plasmatic {

Line::Line(const std::array<Integer, 2> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

Float Line::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;

    const auto length = std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2));

    const auto xi = std::sqrt(std::pow(coord.x - x0, 2) + std::pow(coord.y - y0, 2)) / length;

    if (index == 0) {
        return 1.0 - xi;
    }

    if (index == 1) {
        return xi;
    }

    Abort("Invalid shape function index: {}", index);
}

} // namespace plasmatic
