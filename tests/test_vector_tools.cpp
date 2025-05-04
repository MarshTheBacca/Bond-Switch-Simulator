#include "vector_tools.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

class WrapArrayTest
    : public ::testing::TestWithParam<
          std::tuple<std::array<double, 2>, std::array<double, 2>,
                     std::array<double, 2>>> {};

TEST_P(WrapArrayTest, WrapsCorrectly) {
  auto [input, dimensions, expected] = GetParam();
  wrapArray(input, dimensions);
  EXPECT_NEAR(input[0], expected[0], 1e-9);
  EXPECT_NEAR(input[1], expected[1], 1e-9);
}

INSTANTIATE_TEST_SUITE_P(WrapArrayTests, WrapArrayTest,
                         ::testing::Values(
                             // Test case: Within bounds
                             std::make_tuple(std::array<double, 2>{3.5, 4.5},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{3.5, 4.5}),
                             // Test case: Bottom left corner
                             std::make_tuple(std::array<double, 2>{0.0, 0.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{0.0, 0.0}),
                             // Test case: Bottom right corner
                             std::make_tuple(std::array<double, 2>{10.0, 0.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{0.0, 0.0}),
                             // Test case: Top left corner
                             std::make_tuple(std::array<double, 2>{0.0, 10.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{0.0, 0.0}),
                             // Test case: Top right corner
                             std::make_tuple(std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{0.0, 0.0}),
                             // Test case: Out-of-bounds positive
                             std::make_tuple(std::array<double, 2>{12.5, 14.5},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{2.5, 4.5}),
                             // Test case: Out-of-bounds negative
                             std::make_tuple(std::array<double, 2>{-2.5, -5.5},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{7.5, 4.5}),
                             // Test case: Multiple wraps
                             std::make_tuple(std::array<double, 2>{-25.3,
                                                                   156.2},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{4.7, 6.2})));

class PBCArrayTest
    : public ::testing::TestWithParam<
          std::tuple<std::array<double, 2>, std::array<double, 2>,
                     std::array<double, 2>, std::array<double, 2>>> {};

TEST_P(PBCArrayTest, PBCCorrectly) {
  auto [vector1, vector2, dimensions, expected] = GetParam();
  std::array<double, 2> result = pbcArray(vector1, vector2, dimensions);
  EXPECT_NEAR(result[0], expected[0], 1e-9);
  EXPECT_NEAR(result[1], expected[1], 1e-9);
}

INSTANTIATE_TEST_SUITE_P(PBCTests, PBCArrayTest,
                         ::testing::Values(
                             // Test case with no wrapping required
                             std::make_tuple(std::array<double, 2>{1.0, 1.0},
                                             std::array<double, 2>{2.0, 2.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{1.0, 1.0}),
                             // Over half dimension range
                             std::make_tuple(std::array<double, 2>{1.0, 1.0},
                                             std::array<double, 2>{7.0, 2.0},
                                             std::array<double, 2>{10.0, 10.0},
                                             std::array<double, 2>{-4.0,
                                                                   1.0})));

class AverageArrayTest
    : public ::testing::TestWithParam<std::tuple<
          std::vector<std::array<double, 2>>, std::array<double, 2>>> {};

TEST_P(AverageArrayTest, AverageArrayCorrectly) {
  auto [arrays, expected] = GetParam();
  std::array<double, 2> result = arrayAverage(arrays);
  EXPECT_NEAR(result[0], expected[0], 1e-9);
  EXPECT_NEAR(result[1], expected[1], 1e-9);
}

INSTANTIATE_TEST_SUITE_P(
    AverageArrayTests, AverageArrayTest,
    ::testing::Values(
        // Test case with no wrapping required
        std::make_tuple(std::vector<std::array<double, 2>>{{1.0, 1.0},
                                                           {2.0, 2.0}},
                        std::array<double, 2>{1.5, 1.5}),

        std::make_tuple(std::vector<std::array<double, 2>>{{1.0, 1.0},
                                                           {2.0, 2.0},
                                                           {3.0, 3.0}},
                        std::array<double, 2>{2.0, 2.0}),
        // Test case with 10 random arrays
        std::make_tuple(std::vector<std::array<double, 2>>{{1924.168, 1235.692},
                                                           {-12.5, 31.6},
                                                           {-120.98, 0.12},
                                                           {897.25, 175.7},
                                                           {-76.4, -0.6},
                                                           {5.67, 8.19},
                                                           {0.0, 13.072},
                                                           {616.7829, -96.2},
                                                           {-79.29, -825.2},
                                                           {1.8, 8.3}},

                        std::array<double, 2>{315.65009, 55.0674})));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
