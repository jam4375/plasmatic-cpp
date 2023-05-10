#include "interface/ProblemTypes/Mechanical.h"

#include "LinearAlgebra/LinearAlgebra.h"

#include <Eigen/Dense>

#include <iostream>

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
Mechanical::Mechanical(const Input &input) : _input(input), _mesh(input.mesh_filename) {}

void Mechanical::Solve() {
    constexpr auto dimension = 3;

    _mesh.AddVectorField("displacement");

    // Create global stiffness matrix and forcing vector
    Matrix stiffness(3 * _mesh.GetNumNodes(), 3 * _mesh.GetNumNodes());
    Vector forcing(3 * _mesh.GetNumNodes());
    Vector displacement_vec_bcs(3 * _mesh.GetNumNodes());

    auto E = _input.youngs_modulus;
    auto v = _input.poisson_ratio;
    auto constant = E / ((1.0 + v) * (1.0 - 2.0 * v));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
    D(0, 0) = constant * (1.0 - v);
    D(1, 1) = constant * (1.0 - v);
    D(2, 2) = constant * (1.0 - v);
    D(3, 3) = constant * 0.5 * (1.0 - 2.0 * v);
    D(4, 4) = constant * 0.5 * (1.0 - 2.0 * v);
    D(5, 5) = constant * 0.5 * (1.0 - 2.0 * v);
    D(0, 1) = constant * v;
    D(0, 2) = constant * v;
    D(1, 0) = constant * v;
    D(1, 2) = constant * v;
    D(2, 0) = constant * v;
    D(2, 1) = constant * v;

    // Loop over elements and add elemental stiffness matrix and forcing vector into the global ones
    for (Integer element_id = 0; element_id < _mesh.GetNumElements(dimension); ++element_id) {
        auto element = _mesh.GetElement(dimension, element_id);

        for (Integer ii = 0; ii < element->NumNodes(); ++ii) {
            auto row = element->GetNodeIndex(ii);
            for (Integer jj = 0; jj < element->NumNodes(); ++jj) {
                auto col = element->GetNodeIndex(jj);

                auto value = element->Integrate(
                    [element, ii, jj, D](const Coord &pos) -> Eigen::MatrixXd {
                        Eigen::MatrixXd Ba = Eigen::MatrixXd::Zero(6, 3);
                        Ba(0, 0) = element->ShapeFnDerivative(ii, 0, pos);
                        Ba(1, 1) = element->ShapeFnDerivative(ii, 1, pos);
                        Ba(2, 2) = element->ShapeFnDerivative(ii, 2, pos);
                        Ba(3, 0) = element->ShapeFnDerivative(ii, 1, pos);
                        Ba(3, 1) = element->ShapeFnDerivative(ii, 0, pos);
                        Ba(4, 1) = element->ShapeFnDerivative(ii, 2, pos);
                        Ba(4, 2) = element->ShapeFnDerivative(ii, 1, pos);
                        Ba(5, 0) = element->ShapeFnDerivative(ii, 2, pos);
                        Ba(5, 2) = element->ShapeFnDerivative(ii, 0, pos);

                        Eigen::MatrixXd Bb = Eigen::MatrixXd::Zero(6, 3);
                        Bb(0, 0) = element->ShapeFnDerivative(jj, 0, pos);
                        Bb(1, 1) = element->ShapeFnDerivative(jj, 1, pos);
                        Bb(2, 2) = element->ShapeFnDerivative(jj, 2, pos);
                        Bb(3, 0) = element->ShapeFnDerivative(jj, 1, pos);
                        Bb(3, 1) = element->ShapeFnDerivative(jj, 0, pos);
                        Bb(4, 1) = element->ShapeFnDerivative(jj, 2, pos);
                        Bb(4, 2) = element->ShapeFnDerivative(jj, 1, pos);
                        Bb(5, 0) = element->ShapeFnDerivative(jj, 2, pos);
                        Bb(5, 2) = element->ShapeFnDerivative(jj, 0, pos);

                        Eigen::MatrixXd inner_mat = Ba.transpose() * D * Bb;

                        return inner_mat;
                    },
                    3, 3);

                for (Integer kk = 0; kk < 3; ++kk) {
                    for (Integer mm = 0; mm < 3; ++mm) {
                        stiffness.AddValue(3 * row + kk, 3 * col + mm, value(kk, mm));
                    }
                }
            }
        }
    }
    stiffness.Assemble();

    // Set boundary conditions
    constexpr auto bc_dimension = 2;

    for (const auto &[physical_name, bc_value] : _input.dirichlet_bcs) {
        auto element_entities1 = _mesh.GetPhysicalEntity(physical_name, bc_dimension);
        for (const auto &element_entity : element_entities1) {
            auto element_inds = _mesh.GetEntity(bc_dimension, element_entity);
            for (const auto &element_ind : element_inds) {
                auto element = _mesh.GetElement(bc_dimension, element_ind);
                for (Integer ii = 0; ii < element->NumNodes(); ++ii) {
                    auto node_ind = element->GetNodeIndex(ii);
                    for (Integer jj = 0; jj < 3; ++jj) {
                        displacement_vec_bcs.SetValue(3 * node_ind + jj, bc_value[static_cast<size_t>(jj)]);

                        stiffness.SetDirichletBC(3 * node_ind + jj, displacement_vec_bcs, forcing);
                    }
                }
            }
        }
    }

    for (const auto &[physical_name, bc_value] : _input.neumann_bcs) {
        auto element_entities2 = _mesh.GetPhysicalEntity(physical_name, bc_dimension);
        for (const auto &element_entity : element_entities2) {
            auto element_inds = _mesh.GetEntity(bc_dimension, element_entity);
            for (const auto &element_ind : element_inds) {
                auto element = _mesh.GetElement(bc_dimension, element_ind);
                for (Integer ii = 0; ii < element->NumNodes(); ++ii) {
                    auto row = element->GetNodeIndex(ii);

                    auto bc_value_copy = bc_value;
                    for (Integer jj = 0; jj < 3; ++jj) {
                        auto value = element->Integrate([element, ii, bc_value_copy, jj](const Coord &pos) -> Float {
                            return bc_value_copy[static_cast<size_t>(jj)] * element->ShapeFn(ii, pos);
                        });

                        forcing.AddValue(3 * row + jj, value);
                    }
                }
            }
        }
    }
    stiffness.Assemble();
    forcing.Assemble();

    // Solve stiffness matrix/forcing vector equation for displacement
    auto displacement_vec = stiffness.Solve(forcing);

    // Transfer solution to mesh field
    for (Integer ii = 0; ii < displacement_vec.Size() / 3; ++ii) {
        _mesh.VectorFieldSetValue("displacement", ii,
                                  {displacement_vec.GetValue(3 * ii), displacement_vec.GetValue(3 * ii + 1),
                                   displacement_vec.GetValue(3 * ii + 2)});
    }
}

} // namespace plasmatic
