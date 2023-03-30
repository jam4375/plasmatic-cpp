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
}
} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    PetscErrorCode ierr = PetscInitialize(&argc, &argv, "", "");
    if (ierr != 0) {
        std::abort();
    }

    return RUN_ALL_TESTS();
}
