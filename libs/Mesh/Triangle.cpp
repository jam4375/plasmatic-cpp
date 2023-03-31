#include "interface/Mesh/Triangle.h"

namespace plasmatic {

Triangle::Triangle(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

} // namespace plasmatic
