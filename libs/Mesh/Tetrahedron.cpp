#include "interface/Mesh/Tetrahedron.h"

namespace plasmatic {

Tetrahedron::Tetrahedron(const std::array<Integer, 4> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

Float Tetrahedron::ComputeVolume(const Coord &p0, const Coord &p1, const Coord &p2, const Coord &p3) {
    auto volume_6x =
        (p1.x * (p2.y * p3.z - p2.z * p3.y) - p1.y * (p2.x * p3.z - p2.z * p3.x) + p1.z * (p2.x * p3.y - p2.y * p3.x)) -
        p0.x * ((p2.y * p3.z - p2.z * p3.y) - (p1.y * p3.z - p3.y * p1.z) + (p1.y * p2.z - p2.y * p1.z)) +
        p0.y * ((p2.x * p3.z - p2.z * p3.x) - (p1.x * p3.z - p3.x * p1.z) + (p1.x * p2.z - p2.x * p1.z)) -
        p0.z * ((p2.x * p3.y - p3.x * p2.y) - (p1.x * p3.y - p1.y * p3.x) + (p1.x * p2.y - p2.x * p1.y));

    return volume_6x / 6.0;
}

Float Tetrahedron::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto volume = Tetrahedron::ComputeVolume(
        (*_nodes)[static_cast<size_t>(_nodeIndices[0])], (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
        (*_nodes)[static_cast<size_t>(_nodeIndices[2])], (*_nodes)[static_cast<size_t>(_nodeIndices[3])]);

    if (index == 0) {
        return Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[2])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[1])], coord) /
               volume;
    }

    if (index == 1) {
        return Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[2])], coord) /
               volume;
    }

    if (index == 2) {
        return Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[0])], coord) /
               volume;
    }

    if (index == 3) {
        return Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                          (*_nodes)[static_cast<size_t>(_nodeIndices[2])], coord) /
               volume;
    }

    Abort("Invalid shape function index: {}", index);
}

Float Tetrahedron::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] const Coord &coord) const {
    auto x0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].x;
    auto y0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].y;
    auto z0 = (*_nodes)[static_cast<size_t>(_nodeIndices[0])].z;
    auto x1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].x;
    auto y1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].y;
    auto z1 = (*_nodes)[static_cast<size_t>(_nodeIndices[1])].z;
    auto x2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].x;
    auto y2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].y;
    auto z2 = (*_nodes)[static_cast<size_t>(_nodeIndices[2])].z;
    auto x3 = (*_nodes)[static_cast<size_t>(_nodeIndices[3])].x;
    auto y3 = (*_nodes)[static_cast<size_t>(_nodeIndices[3])].y;
    auto z3 = (*_nodes)[static_cast<size_t>(_nodeIndices[3])].z;

    auto volume_6x = (x1 * (y2 * z3 - z2 * y3) - y1 * (x2 * z3 - z2 * x3) + z1 * (x2 * y3 - y2 * x3)) -
                     x0 * ((y2 * z3 - z2 * y3) - (y1 * z3 - y3 * z1) + (y1 * z2 - y2 * z1)) +
                     y0 * ((x2 * z3 - z2 * x3) - (x1 * z3 - x3 * z1) + (x1 * z2 - x2 * z1)) -
                     z0 * ((x2 * y3 - x3 * y2) - (x1 * y3 - y1 * x3) + (x1 * y2 - x2 * y1));

    if (index == 0) {
        auto b0 = y1 * (z2 - z3) - z1 * (y2 - y3) + (y2 * z3 - z2 * y3);
        auto c0 = z1 * (x3 - x2) - (z2 * x3 - x2 * z3) + x1 * (z2 - z3);
        auto d0 = (x2 * y3 - y2 * x3) - x1 * (y3 - y2) + y1 * (x3 - x2);

        if (dimension == 0) {
            return b0 / volume_6x;
        }
        if (dimension == 1) {
            return c0 / volume_6x;
        }
        if (dimension == 2) {
            return d0 / volume_6x;
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    if (index == 1) {
        auto b1 = y2 * (z3 - z0) - z2 * (y3 - y0) + (y3 * z0 - z3 * y0);
        auto c1 = z2 * (x0 - x3) - (z3 * x0 - x3 * z0) + x2 * (z3 - z0);
        auto d1 = (x3 * y0 - y3 * x0) - x2 * (y0 - y3) + y2 * (x0 - x3);

        if (dimension == 0) {
            return b1 / volume_6x;
        }
        if (dimension == 1) {
            return c1 / volume_6x;
        }
        if (dimension == 2) {
            return d1 / volume_6x;
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    if (index == 2) {
        auto b2 = y3 * (z0 - z1) - z3 * (y0 - y1) + (y0 * z1 - z0 * y1);
        auto c2 = z3 * (x1 - x0) - (z0 * x1 - x0 * z1) + x3 * (z0 - z1);
        auto d2 = (x0 * y1 - y0 * x1) - x3 * (y1 - y0) + y3 * (x1 - x0);

        if (dimension == 0) {
            return b2 / volume_6x;
        }
        if (dimension == 1) {
            return c2 / volume_6x;
        }
        if (dimension == 2) {
            return d2 / volume_6x;
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    if (index == 3) {
        auto b3 = y0 * (z1 - z2) - z0 * (y1 - y2) + (y1 * z2 - z1 * y2);
        auto c3 = z0 * (x2 - x1) - (z1 * x2 - x1 * z2) + x0 * (z1 - z2);
        auto d3 = (x1 * y2 - y1 * x2) - x0 * (y2 - y1) + y0 * (x2 - x1);

        if (dimension == 0) {
            return b3 / volume_6x;
        }
        if (dimension == 1) {
            return c3 / volume_6x;
        }
        if (dimension == 2) {
            return d3 / volume_6x;
        }
        Abort("Invalid shape function derivative dimension: {}", dimension);
    }

    Abort("Invalid shape function index: {}", index);
}

Float Tetrahedron::Integrate([[maybe_unused]] const std::function<Float(const Coord &)> integrand) const {
    constexpr auto alpha = 0.5854102;
    constexpr auto beta = 0.1381966;

    std::vector<std::array<Float, 4>> gauss_coords = {
        {alpha, beta, beta, beta}, {beta, alpha, beta, beta}, {beta, beta, alpha, beta}, {beta, beta, beta, alpha}};

    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    std::vector<Float> weights = {0.25, 0.25, 0.25, 0.25};

    Float result = 0.0;

    auto volume = Tetrahedron::ComputeVolume(
        (*_nodes)[static_cast<size_t>(_nodeIndices[0])], (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
        (*_nodes)[static_cast<size_t>(_nodeIndices[2])], (*_nodes)[static_cast<size_t>(_nodeIndices[3])]);

    for (size_t ii = 0; ii < gauss_coords.size(); ++ii) {
        // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
        Coord point = {.x = 0.0, .y = 0.0, .z = 0.0};

        for (size_t jj = 0; jj < 4; ++jj) {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            point.x += gauss_coords[ii][jj] * (*_nodes)[static_cast<size_t>(_nodeIndices[jj])].x;
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            point.y += gauss_coords[ii][jj] * (*_nodes)[static_cast<size_t>(_nodeIndices[jj])].y;
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            point.z += gauss_coords[ii][jj] * (*_nodes)[static_cast<size_t>(_nodeIndices[jj])].z;
        }

        result += weights[ii] * integrand(point) * volume;
    }

    return result;
}

} // namespace plasmatic
