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

TEST(MeshTest, Tetrahedron) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 1.0});

    Tetrahedron tet({0, 1, 2, 3}, nodes);

    for (Integer ii = 0; ii < 4; ++ii) {
        EXPECT_DOUBLE_EQ(tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]), 1.0) << ", ii = " << ii;

        for (Integer jj = 0; jj < 4; ++jj) {
            if (ii == jj) {
                continue;
            }

            EXPECT_DOUBLE_EQ(tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(jj)]), 0.0)
                << ", ii = " << ii << ", jj = " << jj;
        }
    }

    for (Integer ii = 0; ii < 4; ++ii) {
        auto sum = 0.0;
        sum += tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]);

        for (Integer jj = 0; jj < 4; ++jj) {
            if (ii == jj) {
                continue;
            }

            sum += tet.ShapeFn(jj, (*nodes)[static_cast<size_t>(ii)]);
        }

        EXPECT_DOUBLE_EQ(sum, 1.0) << ", ii = " << ii;
    }

    EXPECT_DOUBLE_EQ(tet.Integrate([](const Coord &) -> Float { return 1.0; }), 1.0 / 6.0);
}

TEST(MeshTest, TetrahedronShapeFnDerivatives) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 1.0});

    for (Integer ii = 0; ii < 4; ++ii) {
        Tetrahedron tet({(0 + ii) % 4, (1 + ii) % 4, (2 + ii) % 4, (3 + ii) % 4}, nodes);

        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((4 - ii) % 4, 0, {.x = 0.0, .y = 0.0, .z = 0.0}), -1.0) << "ii = " << ii;
        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((4 - ii) % 4, 1, {.x = 0.0, .y = 0.0, .z = 0.0}), -1.0) << "ii = " << ii;
        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((4 - ii) % 4, 2, {.x = 0.0, .y = 0.0, .z = 0.0}), -1.0) << "ii = " << ii;

        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((5 - ii) % 4, 0, {.x = 1.0, .y = 0.0, .z = 0.0}), 1.0) << "ii = " << ii;
        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((6 - ii) % 4, 1, {.x = 0.0, .y = 1.0, .z = 0.0}), 1.0) << "ii = " << ii;
        EXPECT_DOUBLE_EQ(tet.ShapeFnDerivative((7 - ii) % 4, 2, {.x = 0.0, .y = 0.0, .z = 1.0}), 1.0) << "ii = " << ii;
    }
}

} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
