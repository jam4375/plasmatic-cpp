#include "interface/ProblemTypes/HeatEq3D.h"

#include "LinearAlgebra/LinearAlgebra.h"

#include <iostream>

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
HeatEq3D::HeatEq3D(const Input &input) : _input(input), _mesh(input.mesh_filename) {}

void HeatEq3D::Solve() {
    constexpr auto dimension = 3;

    _mesh.AddScalarField("temperature");

    // Create global stiffness matrix and forcing vector
    Matrix stiffness(_mesh.GetNumNodes(), _mesh.GetNumNodes());
    Vector forcing(_mesh.GetNumNodes());
    Vector temperature_vec_bcs(_mesh.GetNumNodes());

    // Loop over elements and add elemental stiffness matrix and forcing vector into the global ones
    for (Integer element_id = 0; element_id < _mesh.GetNumElements(dimension); ++element_id) {
        auto element = _mesh.GetElement(dimension, element_id);

        for (Integer ii = 0; ii < element->NumNodes(); ++ii) {
            auto row = element->GetNodeIndex(ii);
            for (Integer jj = 0; jj < element->NumNodes(); ++jj) {
                auto col = element->GetNodeIndex(jj);

                auto value = element->Integrate([element, ii, jj, this](const Coord &pos) -> Float {
                    return _input.thermal_conductivity *
                           (element->ShapeFnDerivative(ii, 0, pos) * element->ShapeFnDerivative(jj, 0, pos) +
                            element->ShapeFnDerivative(ii, 1, pos) * element->ShapeFnDerivative(jj, 1, pos) +
                            element->ShapeFnDerivative(ii, 2, pos) * element->ShapeFnDerivative(jj, 2, pos));
                });

                stiffness.AddValue(row, col, value);
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
                    temperature_vec_bcs.SetValue(node_ind, bc_value);

                    stiffness.SetDirichletBC(node_ind, temperature_vec_bcs, forcing);
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
                    auto value = element->Integrate([element, ii, bc_value_copy](const Coord &pos) -> Float {
                        return bc_value_copy * element->ShapeFn(ii, pos);
                    });

                    forcing.AddValue(row, value);
                }
            }
        }
    }
    stiffness.Assemble();
    forcing.Assemble();

    // Solve stiffness matrix/forcing vector equation for temperature
    Log::Info("Beginning linear solve");
    auto temperature_vec = stiffness.Solve(forcing);
    Log::Info("Finished linear solve");

    // Transfer solution to mesh field
    for (Integer ii = 0; ii < temperature_vec.Size(); ++ii) {
        _mesh.ScalarFieldSetValue("temperature", ii, temperature_vec.GetValue(ii));
    }
}

} // namespace plasmatic
