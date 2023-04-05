#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(MeshTest, Simple) {
    auto filename = GetExecutablePath() / "assets/Mesh/mesh2d.msh";
    Mesh mesh(filename);

    EXPECT_EQ(mesh.GetNumNodes(), 11);
    EXPECT_EQ(mesh.GetNumElements(1), 8);
    EXPECT_EQ(mesh.GetNumElements(2), 12);

    mesh.WriteVTK("mesh2d.vtk");
}
} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
