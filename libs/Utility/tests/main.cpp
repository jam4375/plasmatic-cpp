#include "Utility/Utility.h"

#include <gtest/gtest.h>

namespace plasmatic {
TEST(UtilityTest, Types) {
    EXPECT_EQ(sizeof(Float), 8);
    EXPECT_EQ(sizeof(Integer), 4);
}
} // namespace plasmatic

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
