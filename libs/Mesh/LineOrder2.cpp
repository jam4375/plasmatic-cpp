#include "interface/Mesh/LineOrder2.h"

#include <Eigen/Dense>

namespace plasmatic {

LineOrder2::LineOrder2(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

Float LineOrder2::ComputeLength(const Coord &p0, const Coord &p1) {
    return std::sqrt(std::pow(p1.x - p0.x, 2) + std::pow(p1.y - p0.y, 2) + std::pow(p1.z - p0.z, 2));
}

Float LineOrder2::PhysicalToParentCoords(const Coord &coord) const {
    const auto xi =
        ComputeLength((*_nodes)[static_cast<size_t>(_nodeIndices[0])], coord) /
        ComputeLength((*_nodes)[static_cast<size_t>(_nodeIndices[0])], (*_nodes)[static_cast<size_t>(_nodeIndices[1])]);
    return xi;
}

Coord LineOrder2::ParentToPhysicalCoords(const Float &xi) const {
    std::array<Float, 2> lambda = {0.0};

    lambda[0] = 1.0 - xi;
    lambda[1] = xi;

    // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
    Coord point = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (size_t ii = 0; ii < lambda.size(); ++ii) {
        point.x += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x * lambda[ii];
        point.y += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y * lambda[ii];
        point.z += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].z * lambda[ii];
    }

    return point;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float LineOrder2::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] Float xi) const {
    Check(dimension == 0, "LineOrder2: Invalid shape function derivative dimension = {}", dimension);

    // lambda_i are the length coordinates, xi are the parent coordinates
    const auto lambda_1 = 1.0 - xi;
    const auto lambda_2 = xi;

    // N_a = (2*lambda_a - 1)*lambda_a     for    a = 1,2
    // N_3 = 4*lambda_1*lambda_2

    if (index == 0) {
        return -1.0 * (4.0 * lambda_1 - 1.0);
    }

    if (index == 1) {
        return 4.0 * lambda_2 - 1.0;
    }

    if (index == 2) {
        // N_3 = 4*lambda_1*lambda_2
        return 4.0 * lambda_2 * -1.0 + 4.0 * lambda_1;
    }

    Abort("Invalid index in ShapeFnDerivative");
}

Float LineOrder2::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto xi = PhysicalToParentCoords(coord);

    const auto lambda_1 = 1.0 - xi;
    const auto lambda_2 = xi;

    if (index == 0) {
        return (2.0 * lambda_1 - 1) * lambda_1;
    }

    if (index == 1) {
        return (2.0 * lambda_2 - 1) * lambda_2;
    }

    if (index == 2) {
        return 4.0 * lambda_1 * lambda_2;
    }

    Abort("Invalid shape function index: {}", index);
}

Float LineOrder2::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] const Coord &coord) const {
    const auto xi = PhysicalToParentCoords(coord);

    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(1, 1);

    for (size_t ii = 0; ii < static_cast<size_t>(this->NumNodes()); ++ii) {
        jacobian(0, 0) +=
            (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x * ShapeFnDerivative(static_cast<Integer>(ii), 0, xi);
    }

    Eigen::VectorXd shape_fn_derivs = Eigen::VectorXd::Zero(1);

    shape_fn_derivs(0) = ShapeFnDerivative(index, 0, xi);

    Eigen::VectorXd global_derivs = jacobian.colPivHouseholderQr().solve(shape_fn_derivs);

    Check(dimension == 0, "LineOrder2: Invalid shape function derivative dimension: {}", dimension);

    return global_derivs(dimension);
}

Float LineOrder2::Integrate(const std::function<Float(const Coord &)> integrand) const {
    std::vector<Float> gauss_coords = {1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)};

    std::vector<Float> weights = {0.5, 0.5};

    Float result = 0.0;

    for (size_t ii = 0; ii < gauss_coords.size(); ++ii) {
        auto gauss_point = ParentToPhysicalCoords(gauss_coords[ii]);

        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(1, 1);
        for (size_t kk = 0; kk < static_cast<size_t>(this->NumNodes()); ++kk) {
            for (Integer jj = 0; jj < 1; ++jj) {
                jacobian(jj, 0) += (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].x *
                                   ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[ii]);
            }
        }

        result += weights[ii] * integrand(gauss_point) * std::abs(jacobian.determinant());
    }

    return result;
}

} // namespace plasmatic
