#include "interface/Mesh/Triangle.h"

namespace plasmatic {

Triangle::Triangle(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

Float Triangle::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto x2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].x;
    auto y2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].y;

    auto area_2x_inv = 1.0 / (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));

    auto lambda0 = area_2x_inv * ((x1 * y2 - x2 * y1) + (y1 - y2) * coord.x + (x2 - x1) * coord.y);
    auto lambda1 = area_2x_inv * ((x2 * y0 - x0 * y2) + (y2 - y0) * coord.x + (x0 - x2) * coord.y);
    auto lambda2 = area_2x_inv * ((x0 * y1 - x1 * y0) + (y0 - y1) * coord.x + (x1 - x0) * coord.y);

    if (index == 0) {
        return lambda0;
    }

    if (index == 1) {
        return lambda1;
    }

    if (index == 2) {
        return lambda2;
    }

    Abort("Invalid shape function index: {}", index);
}

Float Triangle::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto x2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].x;
    auto y2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].y;

    auto area_2x_inv = 1.0 / (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));

    if (index == 0) {
        if (dimension == 0) {
            return area_2x_inv * (y1 - y2);
        }
        if (dimension == 1) {
            return area_2x_inv * (x2 - x1);
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    if (index == 1) {
        if (dimension == 0) {
            return area_2x_inv * (y2 - y0);
        }
        if (dimension == 1) {
            return area_2x_inv * (x0 - x2);
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    if (index == 2) {
        if (dimension == 0) {
            return area_2x_inv * (y0 - y1);
        }
        if (dimension == 1) {
            return area_2x_inv * (x1 - x0);
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    Abort("Invalid shape function index: {}", index);
}

Float Triangle::Integrate(const std::function<Float(const Coord &)> integrand) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto x2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].x;
    auto y2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].y;

    // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
    Coord center = {.x = (x0 + x1 + x2) / 3.0, .y = (y0 + y1 + y2) / 3.0};

    auto area = 0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));

    constexpr auto weight = 1.0;

    return weight * area * integrand(center);
}

} // namespace plasmatic
