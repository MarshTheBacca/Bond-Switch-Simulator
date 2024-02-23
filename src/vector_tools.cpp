#include "vector_tools.h"
/**
 * @brief Calculates the difference between two vectors with periodic boundary conditions
 * @param vector1 First vector
 * @param vector2 Second vector
 * @param dimensions Dimensions of the system xhi, yhi, (zhi if 3D)
 * @throw std::invalid_argument if the sizes of the vectors do not match the dimensions
 * @return The difference vector
 */
std::vector<double> pbcVector(const std::vector<double> &vector1, const std::vector<double> &vector2, const std::vector<double> &dimensions) {
    if (vector1.size() != vector2.size() || vector1.size() != dimensions.size()) {
        std::ostringstream oss;
        oss << "Vector sizes do not match in PBC function: " << vector1.size() << " " << vector2.size() << " " << dimensions.size();
        throw std::invalid_argument(oss.str());
    }
    std::vector<double> differenceVector;
    differenceVector.reserve(vector1.size());

    for (size_t i = 0; i < vector1.size(); ++i) {
        double dimension_range = dimensions[i];
        double halfDimensionRange = dimension_range / 2;
        double difference = vector2[i] - vector1[i];
        if (difference > halfDimensionRange) {
            difference -= dimension_range;
        } else if (difference < -halfDimensionRange) {
            difference += dimension_range;
        }
        differenceVector.emplace_back(difference);
    }
    return differenceVector;
}