#include "vector_tools.h"
#include <cmath>

/**
 * @brief Normalises the values of a map in place
 * @tparam T The type of the map
 */
void normaliseMap(std::map<size_t, double> &map) {
  double sum =
      std::accumulate(map.begin(), map.end(), 0.0,
                      [](const double &a, const std::pair<size_t, double> &b) {
                        return a + b.second;
                      });
  for (auto &[key, value] : map) {
    value /= sum;
  }
}

/**
 * @brief Prints a nested map to the console
 * @param map The map to print
 */
void showNestedMap(const std::map<int, std::map<int, double>> &map) {
  for (const auto &[key, innerMap] : map) {
    std::cout << key << ": {";
    for (const auto &[innerKey, value] : innerMap) {
      std::cout << innerKey << ": " << value << ", ";
    }
    std::cout << "}" << std::endl;
  }
}

/**
 * @brief Calculates angle between two vectors
 * @param vector1 First vector
 * @param vector2 Second vector
 * @return Angle between the two vectors in radians
 */
double getClockwiseAngleBetweenVectors(const std::array<double, 2> &vector1,
                                       const std::array<double, 2> &vector2) {
  double dotProduct = vector1[0] * vector2[0] + vector1[1] * vector2[1];
  double magnitudeProduct =
      std::sqrt(vector1[0] * vector1[0] + vector1[1] * vector1[1]) *
      std::sqrt(vector2[0] * vector2[0] + vector2[1] * vector2[1]);
  double angle = std::acos(dotProduct / magnitudeProduct);

  // If the cross product is positive, subtract the angle from 2π to get the
  // angle in the range [π, 2π]
  if (vector1[0] * vector2[1] - vector1[1] * vector2[0] > 0) {
    angle = 2 * M_PI - angle;
  }

  return angle;
}

/**
 * @brief Get the clockwise angle of the vector between two nodes relative to
 * the x axis taking into account periodic boundary conditions
 * @param coord1 Coordinates of the first node
 * @param coord2 Coordinates of the second node
 * @param dimensions Dimensions of the periodic box
 * @return Angle of the vector between the two nodes relative to the x axis in
 * radians
 */
double getClockwiseAnglePBC(const std::array<double, 2> &coord1,
                            const std::array<double, 2> &coord2,
                            const std::array<double, 2> &dimensions) {
  std::array<double, 2> relativeCoord = pbcArray(coord1, coord2, dimensions);
  double angle = std::atan2(relativeCoord[1], relativeCoord[0]);
  if (angle < 0) {
    angle += 2 * M_PI;
  }
  return 2 * M_PI - angle;
}

/**
 * @brief Get the clockwise angle of a point relative to the x axis
 * @param point Coordinates of the point
 */
double getClockwiseAngle(const std::array<double, 2> &point) {
  double angle = std::atan2(point[1], point[0]);
  if (angle < 0) {
    angle += 2 * M_PI;
  }
  return 2 * M_PI - angle;
}

/**
 * @brief Sorts a vector of coordinates in clockwise order relative to [0, 0]
 * @param coords Vector of coordinates to sort
 */
void sortCoordinatesClockwise(std::vector<std::array<double, 2>> &coords) {
  std::ranges::sort(coords, [](const std::array<double, 2> &a,
                               const std::array<double, 2> &b) {
    double angleA = getClockwiseAngle(a);
    double angleB = getClockwiseAngle(b);
    return angleA < angleB;
  });
}

/**
 * @brief Calculates the area of a polygon given its vertices
 * @param vertices Vertices of the polygon
 * @return Area of the polygon
 */
double calculatePolygonArea(std::vector<std::array<double, 2>> &vertices) {
  if (vertices.empty()) {
    return 0.0;
  }
  sortCoordinatesClockwise(vertices);
  double area = 0.0;
  for (size_t i = 0; i < vertices.size() - 1; ++i) {
    area += vertices[i][0] * vertices[i + 1][1] -
            vertices[i + 1][0] * vertices[i][1];
  }
  area += vertices.back()[0] * vertices.front()[1] -
          vertices.front()[0] * vertices.back()[1];
  area /= 2.0;
  return std::abs(area);
}