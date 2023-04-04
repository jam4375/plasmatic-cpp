#include "ProblemTypes/ProblemTypes.h"
#include "LinearAlgebra/LinearAlgebra.h"

#include <gtest/gtest.h>

namespace plasmatic {

TEST(ProblemTypesTest, Simple) {
    auto filename = GetExecutablePath() / "assets/Mesh/mesh2d.msh";

    ProblemType problem(filename);

    problem.Solve();

    problem.WriteVTK("heat2d.vtk");
}

} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    PetscErrorCode ierr = PetscInitialize(&argc, &argv, "", "");
    plasmatic::Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return RUN_ALL_TESTS();
}
