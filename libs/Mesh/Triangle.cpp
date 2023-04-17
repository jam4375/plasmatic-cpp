#include "interface/Mesh/Triangle.h"

#include <Eigen/Dense>

namespace plasmatic {

Triangle::Triangle(const std::array<Integer, 3> &node_indices, const std::shared_ptr<std::vector<Coord>> &nodes)
    : _nodeIndices(node_indices), _nodes(nodes) {}

Float Triangle::ComputeArea(const Coord &p1, const Coord &p2, const Coord &p3) {
    auto x12 = p2.x - p1.x;
    auto y12 = p2.y - p1.y;
    auto z12 = p2.z - p1.z;

    auto x13 = p3.x - p1.x;
    auto y13 = p3.y - p1.y;
    auto z13 = p3.z - p1.z;

    return 0.5 * std::sqrt(std::pow(y12 * z13 - z12 * y13, 2) + std::pow(z12 * x13 - x12 * z13, 2) +
                           std::pow(x12 * y13 - y12 * x13, 2));
}

std::array<Float, 2> Triangle::PhysicalToParentCoords(const Coord &coord) const {
    auto area = Triangle::ComputeArea((*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                      (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                      (*_nodes)[static_cast<size_t>(_nodeIndices[2])]);

    auto lambda2 = Triangle::ComputeArea(coord, (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                         (*_nodes)[static_cast<size_t>(_nodeIndices[2])]) /
                   area;
    auto lambda3 = Triangle::ComputeArea(coord, (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                         (*_nodes)[static_cast<size_t>(_nodeIndices[1])]) /
                   area;

    const auto xi = lambda2;
    const auto eta = lambda3;

    return {xi, eta};
}

Coord Triangle::ParentToPhysicalCoords(const std::array<Float, 2> &parent_coords) const {
    std::array<Float, 3> lambda = {0.0};

    const auto &xi = parent_coords[0];
    const auto &eta = parent_coords[1];

    lambda[0] = 1.0 - xi - eta;
    lambda[1] = xi;
    lambda[2] = eta;

    // NOLINTNEXTLINE(clang-diagnostic-pre-c++20-compat-pedantic)
    Coord point = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (size_t ii = 0; ii < static_cast<size_t>(this->NumNodes()); ++ii) {
        point.x += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x * lambda[ii];
        point.y += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y * lambda[ii];
        point.z += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].z * lambda[ii];
    }

    return point;
}

Float Triangle::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    auto parent_coords = PhysicalToParentCoords(coord);

    const auto &xi = parent_coords[0];
    const auto &eta = parent_coords[1];

    if (index == 0) {
        return 1.0 - xi - eta;
    }

    if (index == 1) {
        return xi;
    }

    if (index == 2) {
        return eta;
    }

    Abort("Invalid shape function index: {}", index);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float Triangle::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] Float xi,
                                  [[maybe_unused]] Float eta) const {
    Check(dimension >= 0 && dimension <= 2, "Invalid dimension in ShapeFnDerivative");

    // lambda_i are the shape functions, xi and eta are the parent coordinates
    // lambda_1 = 1 - xi - eta
    // lambda_2 = xi
    // lambda_3 = eta

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

    Abort("Invalid index in ShapeFnDerivative");
}

Float Triangle::ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const {
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(2, 2);

    auto parent_coords = PhysicalToParentCoords(coord);

    for (size_t ii = 0; ii < 3; ++ii) {
        for (Integer jj = 0; jj < 2; ++jj) {
            jacobian(jj, 0) += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x *
                               ShapeFnDerivative(static_cast<Integer>(ii), jj, parent_coords[0], parent_coords[1]);
            jacobian(jj, 1) += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y *
                               ShapeFnDerivative(static_cast<Integer>(ii), jj, parent_coords[0], parent_coords[1]);
        }
    }

    Eigen::VectorXd shape_fn_derivs = Eigen::VectorXd::Zero(2);

    for (Integer jj = 0; jj < 2; ++jj) {
        shape_fn_derivs(jj) = ShapeFnDerivative(index, jj, parent_coords[0], parent_coords[1]);
    }

    Eigen::VectorXd global_derivs = jacobian.colPivHouseholderQr().solve(shape_fn_derivs);

    Check(dimension >= 0 && dimension <= 2, "Invalid shape function derivative dimension: {}", dimension);
    return global_derivs(dimension);
}

Float Triangle::Integrate(const std::function<Float(const Coord &)> integrand) const {
    std::vector<std::array<Float, 2>> gauss_coords = {{1.0 / 3.0, 1.0 / 3.0}};

    std::vector<Float> weights = {1.0};

    Float result = 0.0;

    for (size_t ii = 0; ii < gauss_coords.size(); ++ii) {
        auto gauss_point = ParentToPhysicalCoords(gauss_coords[ii]);

        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(2, 2);
        for (size_t kk = 0; kk < static_cast<size_t>(this->NumNodes()); ++kk) {
            for (Integer jj = 0; jj < 2; ++jj) {
                jacobian(jj, 0) +=
                    (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].x *
                    ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[kk][0], gauss_coords[kk][1]);
                jacobian(jj, 1) +=
                    (*_nodes)[static_cast<size_t>(_nodeIndices[kk])].y *
                    ShapeFnDerivative(static_cast<Integer>(kk), jj, gauss_coords[kk][0], gauss_coords[kk][1]);
            }
        }

        result += weights[ii] * integrand(gauss_point) * (0.5 * jacobian.determinant());
    }

    return result;
}

} // namespace plasmatic
