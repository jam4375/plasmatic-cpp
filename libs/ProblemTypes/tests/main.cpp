#include "LinearAlgebra/LinearAlgebra.h"
#include "ProblemTypes/ProblemTypes.h"

#include <gtest/gtest.h>

namespace plasmatic {

TEST(ProblemTypesTest, HeatEq2D) {
    HeatEq2D::Input input = {.mesh_filename = GetExecutablePath() / "assets/ProblemTypes/mesh2d.msh",
                             .thermal_conductivity = 1.0,
                             .dirichlet_bcs = {{"physical_curve_1", 100.0}},
                             .neumann_bcs = {{"physical_curve_2", -100.0}}};

    HeatEq2D problem(input);

    problem.Solve();

    problem.WriteVTK("heat2d.vtk");
}

TEST(ProblemTypesTest, HeatEq2D_quadratic_mesh) {
    HeatEq2D::Input input = {.mesh_filename = GetExecutablePath() / "assets/ProblemTypes/mesh2d_quadratic.msh",
                             .thermal_conductivity = 1.0,
                             .dirichlet_bcs = {{"physical_curve_1", 100.0}},
                             .neumann_bcs = {{"physical_curve_2", -100.0}}};

    HeatEq2D problem(input);

    problem.Solve();

    problem.WriteVTK("heat2d_quadratic.vtk");
}

TEST(ProblemTypesTest, HeatEq3D) {
    HeatEq3D::Input input = {.mesh_filename = GetExecutablePath() / "assets/ProblemTypes/mesh3d.msh",
                             .thermal_conductivity = 1.0,
                             .dirichlet_bcs = {{"fixed", 100.0}, {"load", -100.0}},
                             .neumann_bcs = {}};

    HeatEq3D problem(input);

    problem.Solve();

    problem.WriteVTK("heat3d.vtk");
}

TEST(ProblemTypesTest, HeatEq3D_quadratic_mesh) {
    HeatEq3D::Input input = {.mesh_filename = GetExecutablePath() / "assets/ProblemTypes/mesh3d_quadratic.msh",
                             .thermal_conductivity = 1.0,
                             .dirichlet_bcs = {{"fixed", 100.0}, {"load", -100.0}},
                             .neumann_bcs = {}};

    HeatEq3D problem(input);

    problem.Solve();

    problem.WriteVTK("heat3d_quadratic.vtk");
}

TEST(ProblemTypesTest, Mechanical) {
    Mechanical::Input input = {.mesh_filename = GetExecutablePath() / "assets/ProblemTypes/mesh3d_quadratic.msh",
                               .youngs_modulus = 69.0e9,
                               .poisson_ratio = 0.32,
                               .dirichlet_bcs = {{"fixed", {0.0, 0.0, 0.0}}},
                               .neumann_bcs = {{"load", {0.0, -100.0, 0.0}}}};

    Mechanical problem(input);

    problem.Solve();

    problem.WriteVTK("mechanical.vtk");
}

} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    PetscErrorCode ierr = PetscInitialize(&argc, &argv, "", "");
    plasmatic::Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return RUN_ALL_TESTS();
}
