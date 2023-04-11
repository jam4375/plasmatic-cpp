#include "interface/Mesh/Line.h"

namespace plasmatic {

Line::Line(const std::array<Integer, 2> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float Line::ShapeFn(Integer index, Float xi) const {
    if (index == 0) {
        return 1.0 - xi;
    }

    if (index == 1) {
        return xi;
    }

    Abort("Invalid shape function index: {}", index);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float Line::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] Float xi) const {
    Check(dimension == 0, "Invalid shape function derivative dimension = {}", dimension);

    if (index == 0) {
        return -1.0;
    }

    if (index == 1) {
        return 1.0;
    }

    Abort("Invalid shape function index: {}", index);
}

Float Line::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto z0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].z;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto z1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].z;

    const auto length = std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2));

    const auto xi =
        std::sqrt(std::pow(coord.x - x0, 2) + std::pow(coord.y - y0, 2) + std::pow(coord.z - z0, 2)) / length;

    return ShapeFn(index, xi);
}

Float Line::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto z0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].z;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto z1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].z;

    const auto length = std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2));

    const auto xi =
        std::sqrt(std::pow(coord.x - x0, 2) + std::pow(coord.y - y0, 2) + std::pow(coord.z - z0, 2)) / length;

    const auto theta = std::atan2(y1 - y0, x1 - x0);

    constexpr auto tol = 1.0e-10;

    if (dimension == 0) {
        if (std::cos(theta) > tol) {
            return 1.0 / (length * std::cos(theta)) * ShapeFnDerivative(index, 0, xi);
        }

        return 0.0;
    }
    if (dimension == 1) {
        if (std::sin(theta) > tol) {
            return 1.0 / (length * std::sin(theta)) * ShapeFnDerivative(index, 0, xi);
        }

        return 0.0;
    }
    Abort("Invalid shape function derivative dimension: {}", dimension);
}

Float Line::Integrate(const std::function<Float(const Coord &)> integrand) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto z0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].z;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto z1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].z;

    // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
    Coord midpoint = {.x = 0.5 * (x0 + x1), .y = 0.5 * (y0 + y1), .z = 0.5 * (z0 + z1)};

    const auto length = std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) + std::pow(z1 - z0, 2));

    constexpr auto weight = 1.0;

    return weight * length * integrand(midpoint);
}

} // namespace plasmatic
