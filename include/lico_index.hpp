// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <climits>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <piecewise_linear_model.hpp>
#include <config.hpp>
#if RESIDUAL_COMPRESS
#include <lico_fastpfor.hpp>
#endif


namespace lico
{
    #define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

    template <typename K, uint64_t Epsilon = 64>
    class LICO {

    public:
        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        uint64_t n;                           ///< The number of elements this index was built on.
        uint32_t Epsilon_Data;                ///< The epsilon value used to build the index.actual it's uint32_t
        uint8_t bpc_epsilon = BIT_WIDTH(Epsilon);

        std::vector<Intercept_Value> seg_intercept;
        std::vector<uint8_t> seg_slope_exponent;
        std::vector<int64_t> seg_slope_significand;
        std::vector<Covered_Value> seg_covered;

        sdsl::int_vector<64> seg_intercept_compress;
        sdsl::int_vector<64> seg_slope_exponent_compress;
        sdsl::int_vector<64> seg_slope_significand_compress;
        sdsl::int_vector<64> seg_covered_compress;

        sdsl::int_vector<64> corrections_none; // corrections for saving, each value <= (epsilon / 4)
        sdsl::bit_vector signs_none; // signs compress for saving
        std::vector<uint32_t> corrections_compress;

        std::vector<Correction_Value> corrections_vector; // corrections for decode

        uint32_t segments_size;
        uint8_t bpc_covered, bpc_intercept, bpc_slope_significand, bpc_slope_exponent;

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();

        using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;

        using canonical_segment = typename internal::OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;

        uint64_t errorPointCount;

        LICO(): n(0), errorPointCount(0) {}

        explicit LICO(const std::vector<K>& data) : LICO(data.begin(), data.end()) {}

        template <typename RandomIt>
        LICO(RandomIt begin, RandomIt end):
            n(std::distance(begin, end)),
            errorPointCount(0){
            Epsilon_Data = Epsilon;
            corrections_vector.resize(n);

#if RESIDUAL_COMPRESS
#else
            corrections_none.resize(((n * bpc_epsilon + 63) >> 6)); // set 3 to avoid insEufficient memory
            signs_none.resize(n);
#endif

            build(begin, end, Epsilon);

#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor")
                corrections_compress = compress_residuals_fastpfor(corrections_vector);
#endif

            segments_compress();

        }

        template <typename RandomIt>
        void build(RandomIt begin, RandomIt end, uint64_t epsilon) {

            auto n = (uint64_t) std::distance(begin, end);
            if (n == 0) return;

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            std::vector<canonical_segment> canonical_segments;
            canonical_segments.reserve(epsilon > 0 ? n / (epsilon * epsilon) : n / 8);

            auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
            auto out_fun = [&canonical_segments](auto cs) { canonical_segments.push_back(cs); };
            auto n_segments = internal::make_segmentation_par(n, epsilon, in_fun, out_fun);
            auto last_n = n_segments;

            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
                auto i = it -> get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it) -> get_first_x() : n;
                build_segments(*it, n, i, j, begin); // build the segment
            }
        }

        template <typename RandomIt>
        void build_segments(const canonical_segment& cs, uint64_t n, uint64_t i, uint64_t j, RandomIt data) {

            uint32_t first = cs.get_first_x();
            if (first == n) return;

            // if (seg_covered.size() == 10) {
            //     std::cerr << "HERE" << std::endl;
            // }

            auto [cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(first, j - i); // fixed point slope and intercept

            if (first == n - 1) {
                seg_intercept.push_back(data[first]);
                seg_slope_exponent.push_back(0);
                seg_slope_significand.push_back(0);
                seg_covered.push_back(1);

#if RESIDUAL_COMPRESS
                corrections_vector[first] = 0;
#else
                set_correction(corrections_none.data(), first, 0, 0, signs_none);
#endif
                return;
            }


            seg_covered.push_back(j - i);
            seg_intercept.push_back(cs_intercept);
            seg_slope_exponent.push_back(cs_exponent);
            seg_slope_significand.push_back(cs_significand);

#if RESIDUAL_COMPRESS
                int64_t last_correction = 0;
                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    int64_t error_diff = error - last_correction;
                    corrections_vector[p] = error_diff;
                    last_correction = error;
                }
#else
                uint64_t* corrections_ptr = corrections_none.data();
                sdsl::bit_vector& signs_ptr = signs_none;

                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    uint8_t sign_value = error >= 0 ? 0 : 1;
                    error = std::abs(error);
                    if (error <= Epsilon)
                        set_correction(corrections_ptr, p, error, sign_value, signs_ptr);
                    else if (error == Epsilon + 1 && sign_value == 0) // use -0 to represent error = epsilon
                        set_correction(corrections_ptr, p, 0, 1, signs_ptr);
                    else
                        std::cerr << "Error: error = -epsilon: " << error << " " << sign_value << std::endl;
                }
#endif
        }

        void segments_compress() {
            segments_size = seg_covered.size();
            auto max_iter_1 = std::max_element(seg_covered.begin(), seg_covered.end());
            bpc_covered = BIT_WIDTH(*max_iter_1);
            auto max_iter_2 = std::max_element(seg_intercept.begin(), seg_intercept.end());
            bpc_intercept = BIT_WIDTH(*max_iter_2);
            auto max_iter_3 = std::max_element(seg_slope_significand.begin(), seg_slope_significand.end());
            bpc_slope_significand = BIT_WIDTH(*max_iter_3);
            auto max_iter_4 = std::max_element(seg_slope_exponent.begin(), seg_slope_exponent.end());
            bpc_slope_exponent = BIT_WIDTH(*max_iter_4);

            seg_covered_compress.resize((segments_size * bpc_covered + 63) >> 6);
            seg_intercept_compress.resize((segments_size * bpc_intercept + 63) >> 6);
            seg_slope_significand_compress.resize((segments_size * bpc_slope_significand + 63) >> 6);
            seg_slope_exponent_compress.resize((segments_size * bpc_slope_exponent + 63) >> 6);

            for (int i = 0; i < seg_covered.size(); i++) {
                auto covered = seg_covered[i];
                set_segment_compress(seg_covered_compress.data(), i, covered, bpc_covered);
                auto intercept = seg_intercept[i];
                set_segment_compress(seg_intercept_compress.data(), i, intercept, bpc_intercept);
                auto slope_significand = seg_slope_significand[i];
                set_segment_compress(seg_slope_significand_compress.data(), i, slope_significand, bpc_slope_significand);
                auto slope_exponent = seg_slope_exponent[i];
                set_segment_compress(seg_slope_exponent_compress.data(), i, slope_exponent, bpc_slope_exponent);
            }
            // clear the original segments
            seg_covered.clear();
            seg_intercept.clear();
            seg_slope_significand.clear();
            seg_slope_exponent.clear();
            corrections_vector.clear();
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Intercept_Value>().swap(seg_intercept);
            std::vector<int64_t>().swap(seg_slope_significand);
            std::vector<uint8_t>().swap(seg_slope_exponent);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        void free_memory() {
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Intercept_Value>().swap(seg_intercept);
            std::vector<int64_t>().swap(seg_slope_significand);
            std::vector<uint8_t>().swap(seg_slope_exponent);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        int64_t seg_approximate(uint32_t i, uint32_t first, uint8_t exponent, int64_t significand, int32_t intercept) const {
            return (int64_t(significand * (i - first)) >> exponent) + intercept;
        }

        uint64_t segment_slope_significand_max() const {
            auto max_iter = std::max_element(seg_slope_significand.begin(), seg_slope_significand.end());
            uint64_t max_slope_significand = *max_iter;
            return max_slope_significand;
        }

        uint32_t segment_slope_exponent_max() const {
            auto max_iter = std::max_element(seg_slope_exponent.begin(), seg_slope_exponent.end());
            uint32_t max_slope_exponent = *max_iter;
            return max_slope_exponent;
        }

        uint64_t size() const {
            return n;
        }
        

        uint64_t total_size_in_bytes() const {
            return segment_size_in_bytes() + corrections_size_in_bytes() + signs_size_in_bytes();
        }

        uint64_t ground_truth_build_size_in_bytes() const {
            return (double(segments_size) * 136.0 / 8.0)  + sizeof(uint8_t) + corrections_size_in_bytes() + signs_size_in_bytes();
        }

        uint64_t segment_size_in_bytes(bool is_la_vector=false) const {
            if (is_la_vector)
                return (double(segments_size) * 128.0 + sizeof(uint8_t)); // for la_vector, one segment need a 64 bit floating slope, 32 bit intercept, 32 bit breakpoint
            else
                return ((seg_covered_compress.bit_size() + seg_intercept_compress.bit_size() + seg_slope_significand_compress.bit_size() + seg_slope_exponent_compress.bit_size()) / CHAR_BIT) + sizeof(uint8_t); //here uint8_t represents the byte need for Epsilon_Data, actually it's only need 8 bit, but we store it as u32 for convenient
        }

        uint64_t corrections_size_in_bytes() const {
#if RESIDUAL_COMPRESS
            return corrections_compress.size() * 4;
#else
            return (corrections_none.bit_size()) / CHAR_BIT;
#endif
        }

        uint64_t signs_size_in_bytes() const {
#if RESIDUAL_COMPRESS
            return 0;
#else
            return (signs_none.bit_size()) / CHAR_BIT;
#endif
        }


        void report_residual_list(std::ofstream &file) {
#if RESIDUAL_COMPRESS
            file << n << std::endl;
            uint32_t first = 0;
            for (uint32_t i = 0; i < segments_size; i++) {
                auto covered = seg_covered[i];
                int32_t last_correction = 0;
                // file << covered << std::endl;
                for (int j = first; j < covered + first; j++) {
                    auto delta_residual = corrections_vector[j];
                    auto true_residual = last_correction + delta_residual;
                    last_correction = true_residual;
                    file << true_residual << " " << delta_residual << " " << zigzag(delta_residual) << std::endl;
                }
                first += covered;
            }
#else
            throw std::runtime_error("don't support to report uncompress residuals");
#endif
        }

        void segment_init() {
            seg_covered.resize(segments_size);
            seg_intercept.resize(segments_size);
            seg_slope_significand.resize(segments_size);
            seg_slope_exponent.resize(segments_size);

            for (int i = 0; i < segments_size; i++) {
                seg_covered[i] = get_segment_compress(seg_covered_compress.data(), i, bpc_covered);
                seg_intercept[i] = get_segment_compress(seg_intercept_compress.data(), i, bpc_intercept);
                seg_slope_significand[i] = get_segment_compress(seg_slope_significand_compress.data(), i, bpc_slope_significand);
                seg_slope_exponent[i] = get_segment_compress(seg_slope_exponent_compress.data(), i, bpc_slope_exponent);
            }
        }

        void residual_init() {
#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor") {
                std::vector<uint32_t> uncompressed_output(n);
                corrections_vector.resize(n);
                decompress_residuals_fastpfor(corrections_compress, n, uncompressed_output, corrections_vector);
            } else {
                throw std::runtime_error("Not support uncompress residuals");
            }
#else
            corrections_vector.resize(n);
            int32_t * corrections_pointer = corrections_vector.data();
            uint32_t first = 0;
            for(uint32_t i = 0; i < segments_size; i++) {
                auto covered = seg_covered[i];
                int32_t last_correction = 0;
                for (Covered_Value j = first; j < first + covered; ++j) {
                    int32_t correction_varify = get_correction(corrections_none.data(), j, signs_none);
                    *corrections_pointer++ = correction_varify - last_correction;
                    last_correction = correction_varify;
                }
                first += covered;
            }
#endif
        }

        void normal_init(){
            segment_init();
            residual_init();
        }

        void normal_clean() {
            seg_covered.clear();
            seg_intercept.clear();
            seg_slope_significand.clear();
            seg_slope_exponent.clear();
            corrections_vector.clear();
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Intercept_Value>().swap(seg_intercept);
            std::vector<int64_t>().swap(seg_slope_significand);
            std::vector<uint8_t>().swap(seg_slope_exponent);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        std::vector<K> normal_decode() {
            std::vector<K> output;
            output.resize(n);
            uint32_t pointer = 0;

            residual_init();

            for (uint32_t i= 0; i < segments_size; i++) {
                auto covered = seg_covered[i];
                auto significand = seg_slope_significand[i];
                auto exponent = seg_slope_exponent[i];
                auto intercept = seg_intercept[i];

                int32_t last_correction = intercept;
                for (int j = 0; j < covered; ++j) {
                    last_correction = last_correction + corrections_vector[pointer];
                    output[pointer] = ((j * significand) >> exponent) + last_correction;
                    pointer++;
                }
            }

            return output;
        }

        uint64_t get_correction_bit_offset(uint64_t i, uint32_t bpc) {
            return i * bpc;
        }

        void set_segment_compress(uint64_t* compress_val, uint64_t i, uint64_t value, uint8_t bpc) {
                uint64_t idx = get_correction_bit_offset(i, bpc);
                sdsl::bits::write_int(compress_val + (idx >> 6u), value, idx & 0x3f, bpc);
            }

        int64_t get_segment_compress(const uint64_t* compress_val, uint64_t i, uint8_t bpc) {
                uint64_t idx = get_correction_bit_offset(i, bpc);
                uint64_t val = sdsl::bits::read_int(compress_val + (idx >> 6u), idx & 0x3f, bpc);
                return val;
            }


        void set_correction(uint64_t* corrections_compress_val, uint64_t i, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_ptr) {
            uint64_t idx = get_correction_bit_offset(i, bpc_epsilon);
            sdsl::bits::write_int(corrections_compress_val + (idx >> 6u), value, idx & 0x3f, bpc_epsilon);
            signs_ptr[i] = sign_value & 0x01;
        }


        int64_t get_correction(const uint64_t* corrections_compress_val, int64_t i, const sdsl::bit_vector& signs_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_epsilon);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_epsilon);
            if (correction == 0 && signs_ptr[i] == 1)
                return Epsilon + 1;
            return signs_ptr[i] == 0 ? correction : -correction;
        }

     };

}
