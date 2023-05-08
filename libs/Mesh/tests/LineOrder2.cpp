#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(MeshTest, LineOrder2) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});

    LineOrder2 line({0, 1, 2}, nodes);

    for (Integer ii = 0; ii < static_cast<Integer>(nodes->size()); ++ii) {
        EXPECT_DOUBLE_EQ(line.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]), 1.0) << ", ii = " << ii;

        for (Integer jj = 0; jj < static_cast<Integer>(nodes->size()); ++jj) {
            if (ii == jj) {
                continue;
            }

            EXPECT_DOUBLE_EQ(line.ShapeFn(ii, (*nodes)[static_cast<size_t>(jj)]), 0.0)
                << ", ii = " << ii << ", jj = " << jj;
        }
    }

    for (Integer ii = 0; ii < static_cast<Integer>(nodes->size()); ++ii) {
        auto sum = 0.0;
        sum += line.ShapeFn(ii, (*nodes)[static_cast<size_t>(ii)]);

        for (Integer jj = 0; jj < static_cast<Integer>(nodes->size()); ++jj) {
            if (ii == jj) {
                continue;
            }

            sum += line.ShapeFn(jj, (*nodes)[static_cast<size_t>(ii)]);
        }

        EXPECT_DOUBLE_EQ(sum, 1.0) << ", ii = " << ii;
    }

    EXPECT_DOUBLE_EQ(line.Integrate([](const Coord &) -> Float { return 1.0; }), 1.0);
}

TEST(MeshTest, LineOrder2ShapeFnDerivatives) {
    auto nodes = std::make_shared<std::vector<Coord>>();

    nodes->push_back({.x = 0.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 1.0, .y = 0.0, .z = 0.0});
    nodes->push_back({.x = 0.5, .y = 0.0, .z = 0.0});

    LineOrder2 line({0, 1, 2}, nodes);

    for (Integer ii = 0; ii < static_cast<Integer>(nodes->size()); ++ii) {
        constexpr auto eps = 1.0e-8;

        Coord p1 = {.x = 0.0, .y = 0.0, .z = 0.0};
        Coord p2x = {.x = p1.x + eps, .y = p1.y, .z = p1.z};

        constexpr auto tol = 1.0e-5;
        {
            auto f1 = line.ShapeFn(ii, p1);
            auto f2 = line.ShapeFn(ii, p2x);

            EXPECT_NEAR(line.ShapeFnDerivative(ii, 0, p1), (f2 - f1) / eps, tol);
        }
    }
}

} // namespace plasmatic
