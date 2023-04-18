#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(MeshTest, TriangleOrder2) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.0});

    TriangleOrder2 tri({0, 1, 2, 3, 4, 5}, nodes);

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        EXPECT_DOUBLE_EQ(tri.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]), 1.0) << ", ii = " << ii;

        for (Integer jj = 0; jj < nodes->size(); ++jj) {
            if (ii == jj) {
                continue;
            }

            EXPECT_DOUBLE_EQ(tri.ShapeFn(ii, (*nodes)[static_cast<size_t>(jj)]), 0.0)
                << ", ii = " << ii << ", jj = " << jj;
        }
    }

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        auto sum = 0.0;
        sum += tri.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]);

        for (Integer jj = 0; jj < nodes->size(); ++jj) {
            if (ii == jj) {
                continue;
            }

            sum += tri.ShapeFn(jj, (*nodes)[static_cast<size_t>(ii)]);
        }

        EXPECT_DOUBLE_EQ(sum, 1.0) << ", ii = " << ii;
    }

    EXPECT_DOUBLE_EQ(tri.Integrate([](const Coord &) -> Float { return 1.0; }), 1.0 / 2.0);
}

TEST(MeshTest, TriangleOrder2ShapeFnDerivatives) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.0});

    TriangleOrder2 tri({0, 1, 2, 3, 4, 5}, nodes);

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        constexpr auto eps = 1.0e-8;

        Coord p1 = {.x = 0.25, .y = 0.1, .z = 0.0};
        Coord p2x = {.x = p1.x + eps, .y = p1.y, .z = p1.z};
        Coord p2y = {.x = p1.x, .y = p1.y + eps, .z = p1.z};

        constexpr auto tol = 1.0e-5;
        {
            auto f1 = tri.ShapeFn(ii, p1);
            auto f2 = tri.ShapeFn(ii, p2x);

            EXPECT_NEAR(tri.ShapeFnDerivative(ii, 0, p1), (f2 - f1) / eps, tol);
        }

        {
            auto f1 = tri.ShapeFn(ii, p1);
            auto f2 = tri.ShapeFn(ii, p2y);

            EXPECT_NEAR(tri.ShapeFnDerivative(ii, 1, p1), (f2 - f1) / eps, tol);
        }
    }
}

} // namespace plasmatic
