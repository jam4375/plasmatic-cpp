#include "ExampleLib/ExampleLib.h"

#include <gtest/gtest.h>

TEST(ExampleLibTest, Simple) {
    EXPECT_TRUE(true);
    EXPECT_NO_THROW(templo::ExampleLib::PrintMessage("this is a test"));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
