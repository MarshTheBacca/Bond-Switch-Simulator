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
  wrapArrayInPlace(input, dimensions);
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
  const auto &[vector1, vector2, dimensions, expected] = GetParam();
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
  const auto &[arrays, expected] = GetParam();
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

class IntersectSetsTest
    : public ::testing::TestWithParam<
          std::tuple<std::unordered_set<int>, std::unordered_set<int>,
                     std::unordered_set<int>>> {};

TEST_P(IntersectSetsTest, IntersectSetsTest) {
  const auto &[set1, set2, expected] = GetParam();
  std::unordered_set<int> result = intersectSets(set1, set2);
  EXPECT_EQ(result, expected);
}

INSTANTIATE_TEST_SUITE_P(IntersectSetsTests, IntersectSetsTest,
                         ::testing::Values(
                             // Test case with no intersection
                             std::make_tuple(std::unordered_set<int>{1, 2, 3},
                                             std::unordered_set<int>{4, 5, 6},
                                             std::unordered_set<int>{}),
                             // Test case with some intersection
                             std::make_tuple(std::unordered_set<int>{1, 2, 3},
                                             std::unordered_set<int>{2, 3, 4},
                                             std::unordered_set<int>{2, 3}),
                             // Test case with full intersection
                             std::make_tuple(std::unordered_set<int>{1, 2, 3},
                                             std::unordered_set<int>{1, 2, 3},
                                             std::unordered_set<int>{1, 2, 3}),
                             std::make_tuple(std::unordered_set<int>{2, 33, 8},
                                             std::unordered_set<int>{8, 33, 2},
                                             std::unordered_set<int>{2, 8,
                                                                     33})));

class CheckAnglesTest
    : public ::testing::TestWithParam<
          std::tuple<std::array<double, 2>, std::vector<std::array<double, 2>>,
                     std::array<double, 2>, double, double, bool>> {};

TEST_P(CheckAnglesTest, CheckAnglesTest) {
  const auto &[coord, connectionCoords, dimensions, minAngle, maxAngle,
               expected] = GetParam();
  bool result =
      checkAnglesPBC(coord, connectionCoords, dimensions, minAngle, maxAngle);
  EXPECT_EQ(result, expected);
}

INSTANTIATE_TEST_SUITE_P(
    CheckAnglesTests, CheckAnglesTest,
    ::testing::Values(
        // Square box with four 90 degree angles
        std::make_tuple(std::array<double, 2>{5.0, 5.0},
                        std::vector<std::array<double, 2>>{
                            {4.0, 5.0}, {6.0, 5.0}, {5.0, 4.0}, {5.0, 6.0}},
                        std::array<double, 2>{10.0, 10.0},
                        (80.0 / 180.0) * M_PI, (100.0 / 180.0) * M_PI, true),

        // An equalateral triangle, should be ~120
        std::make_tuple(std::array<double, 2>{0.0, 0.0},
                        std::vector<std::array<double, 2>>{
                            {1.0, 0.0}, {9.666666, 9.5}, {9.666666, 0.5}},
                        std::array<double, 2>{10.0, 10.0},
                        (100.0 / 180.0) * M_PI, (140.0 / 180.0) * M_PI, true)));

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
