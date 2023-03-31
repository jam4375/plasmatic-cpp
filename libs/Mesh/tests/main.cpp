#include "Mesh/Mesh.h"

#include <gtest/gtest.h>

TEST(MeshTest, Simple) { EXPECT_TRUE(true); }

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
