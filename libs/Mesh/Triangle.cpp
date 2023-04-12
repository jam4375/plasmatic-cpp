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

std::array<Float, 3> Triangle::PhysicalToParentCoords(const Coord &coord) const {
    auto area = Triangle::ComputeArea((*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                      (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                      (*_nodes)[static_cast<size_t>(_nodeIndices[2])]);

    auto lambda1 = Triangle::ComputeArea(coord, (*_nodes)[static_cast<size_t>(_nodeIndices[1])],
                                         (*_nodes)[static_cast<size_t>(_nodeIndices[2])]) /
                   area;
    auto lambda2 = Triangle::ComputeArea(coord, (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                         (*_nodes)[static_cast<size_t>(_nodeIndices[2])]) /
                   area;
    auto lambda3 = Triangle::ComputeArea(coord, (*_nodes)[static_cast<size_t>(_nodeIndices[0])],
                                         (*_nodes)[static_cast<size_t>(_nodeIndices[1])]) /
                   area;

    return {lambda1, lambda2, lambda3};
}

Float Triangle::ShapeFn([[maybe_unused]] Integer index, [[maybe_unused]] const Coord &coord) const {
    Check(index >= 0 && index <= 3, "Invalid shape function index: {}", index);

    return PhysicalToParentCoords(coord)[static_cast<size_t>(index)];
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Float Triangle::ShapeFnDerivative(Integer index, Integer dimension, [[maybe_unused]] Float lambda1,
                                  [[maybe_unused]] Float lambda2) const {
    if (index == 0) {
        if (dimension == 0) {
            return 1.0;
        }

        return 0.0;
    }

    if (index == 1) {
        if (dimension == 1) {
            return 1.0;
        }

        return 0.0;
    }

    if (index == 2) {
        return -1.0;
    }

    Abort("Invalid index in ShapeFnDerivative");
}

Float Triangle::ShapeFnDerivative(Integer index, Integer dimension, const Coord &coord) const {
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(2, 2);

    auto parent_coords = PhysicalToParentCoords(coord);

    for (size_t ii = 0; ii < 3; ++ii) {
        for (Integer jj = 0; jj < 2; ++jj) {
            jacobian(jj, 0) += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].x *
                               ShapeFnDerivative(index, jj, parent_coords[0], parent_coords[1]);
            jacobian(jj, 1) += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].y *
                               ShapeFnDerivative(index, jj, parent_coords[0], parent_coords[1]);
            // jacobian(jj, 2) += (*_nodes)[static_cast<size_t>(_nodeIndices[ii])].z *
            //                     ShapeFnDerivative(index, jj, parent_coords[0], parent_coords[1]);
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
