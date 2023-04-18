#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(MeshTest, TetrahedronOrder2) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 1.0});

    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.5});
    nodes->push_back({.x = 0.5, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.5});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.5});

    TetrahedronOrder2 tet({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, nodes);

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        EXPECT_DOUBLE_EQ(tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]), 1.0) << ", ii = " << ii;

        for (Integer jj = 0; jj < nodes->size(); ++jj) {
            if (ii == jj) {
                continue;
            }

            EXPECT_DOUBLE_EQ(tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(jj)]), 0.0)
                << ", ii = " << ii << ", jj = " << jj;
        }
    }

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        auto sum = 0.0;
        sum += tet.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]);

        for (Integer jj = 0; jj < nodes->size(); ++jj) {
            if (ii == jj) {
                continue;
            }

            sum += tet.ShapeFn(jj, (*nodes)[static_cast<size_t>(ii)]);
        }

        EXPECT_DOUBLE_EQ(sum, 1.0) << ", ii = " << ii;
    }

    EXPECT_DOUBLE_EQ(tet.Integrate([](const Coord &) -> Float { return 1.0; }), 1.0 / 6.0);
}

TEST(MeshTest, TetrahedronOrder2ShapeFnDerivatives) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 1.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 1.0});

    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.5});
    nodes->push_back({.x = 0.5, .y = 0.5, .z = 0.0});
    nodes->push_back({.x = 0.0, .y = 0.5, .z = 0.5});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.5});

    TetrahedronOrder2 tet({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, nodes);

    for (Integer ii = 0; ii < nodes->size(); ++ii) {
        constexpr auto eps = 1.0e-8;

        Coord p1 = {.x = 0.0, .y = 0.0, .z = 0.0};
        Coord p2x = {.x = p1.x + eps, .y = p1.y, .z = p1.z};
        Coord p2y = {.x = p1.x, .y = p1.y + eps, .z = p1.z};
        Coord p2z = {.x = p1.x, .y = p1.y, .z = p1.z + eps};

        constexpr auto tol = 1.0e-5;
        {
            auto f1 = tet.ShapeFn(ii, p1);
            auto f2 = tet.ShapeFn(ii, p2x);

            EXPECT_NEAR(tet.ShapeFnDerivative(ii, 0, p1), (f2 - f1) / eps, tol);
        }
        {
            auto f1 = tet.ShapeFn(ii, p1);
            auto f2 = tet.ShapeFn(ii, p2y);

            EXPECT_NEAR(tet.ShapeFnDerivative(ii, 1, p1), (f2 - f1) / eps, tol);
        }
        {
            auto f1 = tet.ShapeFn(ii, p1);
            auto f2 = tet.ShapeFn(ii, p2z);

            EXPECT_NEAR(tet.ShapeFnDerivative(ii, 2, p1), (f2 - f1) / eps, tol);
        }
    }
}

} // namespace plasmatic
