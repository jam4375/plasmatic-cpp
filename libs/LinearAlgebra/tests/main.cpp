#include "LinearAlgebra/LinearAlgebra.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(LinearAlgebraTest, Vector) {
    Vector vec(5);

    EXPECT_EQ(vec.Size(), 5);

    vec.SetValue(0, 24.0);
    vec.Assemble();
    EXPECT_DOUBLE_EQ(vec.GetValue(0), 24.0);

    vec.AddValue(0, 20.0);
    vec.Assemble();
    EXPECT_DOUBLE_EQ(vec.GetValue(0), 44.0);

    Vector vec2(5);
    vec2.SetValue(0, 2.0);
    vec2.SetValue(1, 8.0);
    vec2.Assemble();

    {
        auto vec3 = vec + vec2;
        EXPECT_DOUBLE_EQ(vec3.GetValue(0), 46.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(1), 8.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(2), 0.0);
    }

    {
        auto vec3 = vec - vec2;
        EXPECT_DOUBLE_EQ(vec3.GetValue(0), 42.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(1), -8.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(2), 0.0);
    }

    {
        auto vec3 = vec;

        vec3 += vec2;
        EXPECT_DOUBLE_EQ(vec3.GetValue(0), 46.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(1), 8.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(2), 0.0);
    }

    {
        auto vec3 = vec;

        vec3 -= vec2;
        EXPECT_DOUBLE_EQ(vec3.GetValue(0), 42.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(1), -8.0);
        EXPECT_DOUBLE_EQ(vec3.GetValue(2), 0.0);
    }
}

TEST(LinearAlgebraTest, Matrix) {
    Matrix mat(5, 5);

    EXPECT_EQ(mat.Rows(), 5);
    EXPECT_EQ(mat.Cols(), 5);

    mat.SetValue(0, 0, 24.0);
    mat.SetValue(1, 1, 0.0);
    mat.Assemble();
    EXPECT_DOUBLE_EQ(mat.GetValue(0, 0), 24.0);

    mat.AddValue(0, 0, 20.0);
    mat.Assemble();
    EXPECT_DOUBLE_EQ(mat.GetValue(0, 0), 44.0);

    Matrix mat2(5, 5);
    mat2.SetValue(0, 0, 2.0);
    mat2.SetValue(1, 1, 8.0);
    mat2.Assemble();

    {
        auto mat3 = mat + mat2;
        mat3.Assemble();
        EXPECT_DOUBLE_EQ(mat3.GetValue(0, 0), 46.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(1, 1), 8.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(2, 0), 0.0);
    }

    {
        auto mat3 = mat - mat2;
        EXPECT_DOUBLE_EQ(mat3.GetValue(0, 0), 42.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(1, 1), -8.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(2, 0), 0.0);
    }

    {
        auto mat3 = mat;

        mat3 += mat2;
        EXPECT_DOUBLE_EQ(mat3.GetValue(0, 0), 46.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(1, 1), 8.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(2, 0), 0.0);
    }

    {
        auto mat3 = mat;

        mat3 -= mat2;
        EXPECT_DOUBLE_EQ(mat3.GetValue(0, 0), 42.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(1, 1), -8.0);
        EXPECT_DOUBLE_EQ(mat3.GetValue(2, 0), 0.0);
    }

    {
        Vector x(5);
        x.SetValue(0, 1);
        x.SetValue(1, 1);
        x.Assemble();

        auto y = mat2 * x;
        EXPECT_DOUBLE_EQ(y.GetValue(0), 2.0);
        EXPECT_DOUBLE_EQ(y.GetValue(1), 8.0);
        EXPECT_DOUBLE_EQ(y.GetValue(2), 0.0);
    }
}

TEST(LinearAlgebraTest, LinearSolver) {
    Matrix mat(5, 5);

    mat.SetValue(0, 0, 1.0);
    mat.SetValue(4, 4, 1.0);
    for (Integer ii = 1; ii < 4; ++ii) {
        mat.SetValue(ii, ii - 1, -1.0);
        mat.SetValue(ii, ii, 2.0);
        mat.SetValue(ii, ii + 1, -1.0);
    }
    mat.Assemble();

    Vector rhs(5);

    for (Integer ii = 0; ii < 5; ++ii) {
        rhs.SetValue(ii, 1.0);
    }

    auto ans = mat.Solve(rhs);

    constexpr auto tol = 1.0e-8;
    EXPECT_NEAR(ans.GetValue(0), 1.0, tol);
    EXPECT_NEAR(ans.GetValue(1), 2.5, tol);
    EXPECT_NEAR(ans.GetValue(2), 3.0, tol);
    EXPECT_NEAR(ans.GetValue(3), 2.5, tol);
    EXPECT_NEAR(ans.GetValue(4), 1.0, tol);
}
} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    PetscErrorCode ierr = PetscInitialize(&argc, &argv, "", "");
    plasmatic::Check(ierr == 0, "PETSc returned a non-zero error code: {}", ierr);

    return RUN_ALL_TESTS();
}
