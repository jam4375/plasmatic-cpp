#include "interface/ProblemTypes/ProblemTypes.h"

#include "LinearAlgebra/LinearAlgebra.h"

#include <iostream>

namespace plasmatic {

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
ProblemType::ProblemType(const std::filesystem::path &mesh_filename) : _mesh(mesh_filename) {}

void ProblemType::Solve() {
    constexpr auto dimension = 2;
    constexpr auto thermal_conductivity = 1.0;

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

                auto value = element->Integrate([element, ii, jj](const Coord &pos) -> Float {
                    return thermal_conductivity *
                           (element->ShapeFnDerivative(ii, 0, pos) * element->ShapeFnDerivative(jj, 0, pos) +
                            element->ShapeFnDerivative(ii, 1, pos) * element->ShapeFnDerivative(jj, 1, pos));
                });

                stiffness.AddValue(row, col, value);
            }
        }
    }
    stiffness.Assemble();

    // Set boundary conditions
    constexpr auto bc1_value = 100.0;
    constexpr auto bc2_value = -100.0;

    auto bc1_node_inds = _mesh.GetNodeEntity(2);
    for (auto node_ind : bc1_node_inds) {
        temperature_vec_bcs.SetValue(node_ind, bc1_value);

        stiffness.SetDirichletBC(node_ind, temperature_vec_bcs, forcing);
    }

    auto bc2_node_inds = _mesh.GetNodeEntity(3);
    for (auto node_ind : bc2_node_inds) {
        temperature_vec_bcs.SetValue(node_ind, bc2_value);

        stiffness.SetDirichletBC(node_ind, temperature_vec_bcs, forcing);
    }
    stiffness.Assemble();
    forcing.Assemble();

    // Solve stiffness matrix/forcing vector equation for temperature
    auto temperature_vec = stiffness.Solve(forcing);

    // Transfer solution to mesh field
    for (Integer ii = 0; ii < temperature_vec.Size(); ++ii) {
        _mesh.ScalarFieldSetValue("temperature", ii, temperature_vec.GetValue(ii));
    }
}

} // namespace plasmatic
