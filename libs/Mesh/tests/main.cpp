#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(MeshTest, Simple) {
    auto filename = GetExecutablePath() / "assets/Mesh/mesh2d.msh";
    Mesh mesh(filename);
}
} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
