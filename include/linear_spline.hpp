#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace ls {

struct Segment {
    int startIndex;
    double slope;           // stored parameter only
};

struct FitResult {
    std::vector<Segment> segments;
    std::size_t segmentBytes;
    std::size_t errorBytes;
    std::int64_t epsilon;
    std::size_t residualViolations; // count of |residual| > epsilon after fitting
};

// Fit a linear spline to a sorted, unique array of int64 values.
// epsilon: maximum absolute error allowed for any point relative to its segment line.
FitResult fitLinearSpline(const std::vector<std::int64_t>& data, std::int64_t epsilon);

struct EncodedSpline {
    std::vector<Segment> segments;          // segment start indices and slopes
    std::int64_t firstValue;                // seed value at index 0
    std::vector<std::int64_t> residuals;    // residuals for all non-start indices in order
    std::size_t n;                          // original array length
    std::int64_t epsilon;                   // max absolute residual bound
};

// Encode: produce segments and residuals (no intercepts, no end indices).
EncodedSpline encodeLinearSpline(const std::vector<std::int64_t>& data, std::int64_t epsilon);

// Decode: reconstruct the original array from EncodedSpline.
std::vector<std::int64_t> decodeLinearSpline(const EncodedSpline& enc);

} // namespace ls
