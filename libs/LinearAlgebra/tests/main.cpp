#include "LinearAlgebra/LinearAlgebra.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(LinearAlgebraTest, Simple) {
    Vector vec(5);

    EXPECT_EQ(vec.Size(), 5);
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
