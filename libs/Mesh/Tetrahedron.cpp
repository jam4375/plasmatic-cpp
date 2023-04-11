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

Float Tetrahedron::ComputeVolumeDerivative(const Coord &p0, const Coord &p1, const Coord &p2, Integer dimension) {
    if (dimension == 0) { // x
        return 1.0 / 6.0 * ((p1.y * p2.z - p1.z * p2.y) + p0.y * (-p2.z + p1.z) - p0.z * (-p2.y + p1.y));
    }

    if (dimension == 1) { // y
        return 1.0 / 6.0 * ((-p1.x * p2.z + p1.z * p2.x) - p0.x * (-p2.z + p1.z) - p0.z * (p2.x - p1.x));
    }

    if (dimension == 2) { // z
        return 1.0 / 6.0 * ((p1.x * p2.y - p1.y * p2.x) - p0.x * (p2.y - p1.y) + p0.y * (p2.x - p1.x));
    }

    Abort("Invalid volume derivative dimension: {}", dimension);
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
    auto volume = Tetrahedron::ComputeVolume(
        (*_nodes)[static_cast<size_t>(_nodeIndices[0])], (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
        (*_nodes)[static_cast<size_t>(_nodeIndices[2])], (*_nodes)[static_cast<size_t>(_nodeIndices[3])]);

    if (index == 0) {
        return Tetrahedron::ComputeVolumeDerivative((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[2])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[1])], dimension) /
               volume;
    }

    if (index == 1) {
        return Tetrahedron::ComputeVolumeDerivative((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[2])], dimension) /
               volume;
    }

    if (index == 2) {
        return Tetrahedron::ComputeVolumeDerivative((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[0])], dimension) /
               volume;
    }

    if (index == 3) {
        return Tetrahedron::ComputeVolumeDerivative((*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                                    (*_nodes)[static_cast<size_t>(_nodeIndices[2])], dimension) /
               volume;
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
