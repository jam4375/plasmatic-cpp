#include "interface/Mesh/Tetrahedron.h"

#include <Eigen/Dense>

namespace plasmatic {

Tetrahedron::Tetrahedron(const std::array<Integer, 4> &node_indices, std::shared_ptr<std::vector<Coord>> nodes)
    : _nodeIndices(node_indices), _nodes(std::move(nodes)) {}

Float Tetrahedron::ComputeVolume(const Coord &p0, const Coord &p1, const Coord &p2, const Coord &p3) {
    auto volume_6x =
        (p1.x * (p2.y * p3.z - p2.z * p3.y) - p1.y * (p2.x * p3.z - p2.z * p3.x) + p1.z * (p2.x * p3.y - p2.y * p3.x)) -
        p0.x * ((p2.y * p3.z - p2.z * p3.y) - (p1.y * p3.z - p3.y * p1.z) + (p1.y * p2.z - p2.y * p1.z)) +
        p0.y * ((p2.x * p3.z - p2.z * p3.x) - (p1.x * p3.z - p3.x * p1.z) + (p1.x * p2.z - p2.x * p1.z)) -
        p0.z * ((p2.x * p3.y - p3.x * p2.y) - (p1.x * p3.y - p1.y * p3.x) + (p1.x * p2.y - p2.x * p1.y));

    return volume_6x / 6.0;
}

std::array<Float, 3> Tetrahedron::PhysicalToParentCoords(const Coord &coord) const {
    auto volume = Tetrahedron::ComputeVolume(
        (*_nodes)[static_cast<size_t>(_nodeIndices[0])], (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
        (*_nodes)[static_cast<size_t>(_nodeIndices[2])], (*_nodes)[static_cast<size_t>(_nodeIndices[3])]);

    auto lambda2 = Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[2])], coord) /
                   volume;

    auto lambda3 = Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[3])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[0])], coord) /
                   volume;

    auto lambda4 = Tetrahedron::ComputeVolume((*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                              (*_nodes)[static_cast<size_t>(_nodeIndices[2])], coord) /
                   volume;

    const auto xi = lambda2;
    const auto eta = lambda3;
    const auto zeta = lambda4;

    return {xi, eta, zeta};
}

Coord Tetrahedron::ParentToPhysicalCoords(const std::array<Float, 3> &parent_coords) const {
    std::array<Float, 4> lambda = {0.0};

    const auto &xi = parent_coords[0];
    const auto &eta = parent_coords[1];
    const auto &zeta = parent_coords[2];

    lambda[0] = 1.0 - xi - eta - zeta;
    lambda[1] = xi;
    lambda[2] = eta;
    lambda[3] = zeta;

    // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
    Coord point = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (size_t ii = 0; ii < static_cast<size_t>(this->NumNodes()); ++ii) {
        point.x += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x * lambda[ii];
        point.y += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y * lambda[ii];
        point.z += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].z * lambda[ii];
    }

    return point;
}

Float Tetrahedron::ShapeFn(Integer index, [[maybe_unused]] const Coord &coord) const {
    auto parent_coords = PhysicalToParentCoords(coord);

    const auto &xi = parent_coords[0];
    const auto &eta = parent_coords[1];
    const auto &zeta = parent_coords[2];

    if (index == 0) {
        return 1.0 - xi - eta - zeta;
    }

    if (index == 1) {
        return xi;
    }

    if (index == 2) {
        return eta;
    }

    if (index == 3) {
        return zeta;
    }

    Abort("Invalid shape function index: {}", index);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float Tetrahedron::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] Float xi,
                                     [[maybe_unused]] Float eta, [[maybe_unused]] Float zeta) const {
    Check(dimension >= 0 && dimension <= 3, "Invalid dimension in ShapeFnDerivative");

    // lambda_i are the shape functions, xi, eta, and zeta are the parent coordinates
    // lambda_1 = 1 - xi - eta - zeta
    // lambda_2 = xi
    // lambda_3 = eta
    // lambda_4 = zeta

    if (index == 0) {
        return -1.0;
    }

    if (index == 1) {
        if (dimension == 0) {
            return 1.0;
        }

        return 0.0;
    }

    if (index == 2) {
        if (dimension == 1) {
            return 1.0;
        }

        return 0.0;
    }

    if (index == 3) {
        if (dimension == 2) {
            return 1.0;
        }

        return 0.0;
    }

    Abort("Invalid index in ShapeFnDerivative");
}

Float Tetrahedron::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] const Coord &coord) const {
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(3, 3);

    auto parent_coords = PhysicalToParentCoords(coord);

    for (size_t ii = 0; ii < 4; ++ii) {
        for (Integer jj = 0; jj < 3; ++jj) {
            jacobian(jj, 0) +=
                (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x *
                ShapeFnDerivative(static_cast<Integer>(ii), jj, parent_coords[0], parent_coords[1], parent_coords[2]);
            jacobian(jj, 1) +=
                (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y *
                ShapeFnDerivative(static_cast<Integer>(ii), jj, parent_coords[0], parent_coords[1], parent_coords[2]);
            jacobian(jj, 2) +=
                (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].z *
                ShapeFnDerivative(static_cast<Integer>(ii), jj, parent_coords[0], parent_coords[1], parent_coords[2]);
        }
    }

    Eigen::VectorXd shape_fn_derivs = Eigen::VectorXd::Zero(3);

    for (Integer jj = 0; jj < 3; ++jj) {
        shape_fn_derivs(jj) = ShapeFnDerivative(index, jj, parent_coords[0], parent_coords[1], parent_coords[2]);
    }

    Eigen::VectorXd global_derivs = jacobian.colPivHouseholderQr().solve(shape_fn_derivs);

    Check(dimension >= 0 && dimension <= 3, "Invalid shape function derivative dimension: {}", dimension);
    return global_derivs(dimension);
}

Float Tetrahedron::Integrate([[maybe_unused]] const std::function<Float(const Coord &)> integrand) const {
    constexpr auto alpha = 0.5854102;
    constexpr auto beta = 0.1381966;

    std::vector<std::array<Float, 3>> gauss_coords = {
        {beta, beta, beta}, {alpha, beta, beta}, {beta, alpha, beta}, {beta, beta, alpha}};

    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    std::vector<Float> weights = {0.25, 0.25, 0.25, 0.25};

    Float result = 0.0;

    for (size_t ii = 0; ii < gauss_coords.size(); ++ii) {
        auto gauss_point = ParentToPhysicalCoords(gauss_coords[ii]);

        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(3, 3);

        for (size_t kk = 0; kk < static_cast<size_t>(this->NumNodes()); ++kk) {
            for (Integer jj = 0; jj < 3; ++jj) {
                jacobian(jj, 0) += (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].x *
                                   ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[ii][0],
                                                     gauss_coords[ii][1], gauss_coords[ii][2]);
                jacobian(jj, 1) += (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].y *
                                   ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[ii][0],
                                                     gauss_coords[ii][1], gauss_coords[ii][2]);
                jacobian(jj, 2) += (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].z *
                                   ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[ii][0],
                                                     gauss_coords[ii][1], gauss_coords[ii][2]);
            }
        }

        result += weights[ii] * integrand(gauss_point) * jacobian.determinant() / 6.0;
    }

    return result;
}

} // namespace plasmatic
