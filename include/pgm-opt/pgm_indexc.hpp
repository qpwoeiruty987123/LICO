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
#include <cmath>
#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <../tools.hpp>
#include <../piecewise_linear_model.hpp>

#include <../rice.hpp>
#include <../fastpfor.hpp>

# ifndef  RESIDUAL_COMPRESS
# define RESIDUAL_COMPRESS 1
# endif

namespace lico
{
    #define BIT_CEIL(x) ((x) < 2 ? 1u : 1u << (64u - __builtin_clzll((x) - 1)))
    #define BIT_WIDTH(x) ((x) == 0 ? 0 : 64 - __builtin_clzll(x))

    template <typename K, uint64_t Epsilon = 64, uint64_t EpsilonRecursive = 0, typename Floating = double>
    class PGM {

    public:

        // static_assert(Epsilon > 0);

        struct Segment;

        typedef int64_t Simd_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        uint64_t n;                           ///< The number of elements this index was built on.
        uint32_t Epsilon_Data;                ///< The epsilon value used to build the index.actual it's uint32_t
        std::vector<uint64_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

        K first_pos;                     ///< The smallest element.
        std::vector<K> seg_first;
        std::vector<Intercept_Value> seg_intercept;
        std::vector<uint8_t> seg_slope_exponent;
        std::vector<int64_t> seg_slope_significand;
        std::vector<Covered_Value> seg_covered;
        sdsl::int_vector<64> seg_first_compress;
        sdsl::int_vector<64> seg_intercept_compress;
        sdsl::int_vector<64> seg_slope_exponent_compress;
        sdsl::int_vector<64> seg_slope_significand_compress;
        sdsl::int_vector<64> seg_covered_compress;


        // orignal delta compress
        uint8_t bit_compress = Epsilon < 7 ? Epsilon < 3 ? 0 : 1 : 2;
        // static constexpr uint8_t bit_compress = 0; // uncompressed
        uint32_t epsilon_compress = Epsilon;
        uint32_t bpc_compress = BIT_WIDTH(Epsilon);
        uint32_t bpc_exception = BIT_WIDTH(Epsilon) + 1; // why + 1 ? because the max of compressed residuals by delta is 2*Epsilon - 1, like 256- (-255)=511. Notably, the max residual of original under Epsilon is Epsilon + 1, we mark it as -0 when it's not second-order compression
        sdsl::bit_vector signs_none; // signs compress for saving
        sdsl::bit_vector signs_exception; // signs exception for saving
        uint64_t corrections_exception_num_write;
        uint64_t corrections_exception_num_read;

        std::vector<uint8_t> corrections_compress_rice; // unused
        std::vector<uint32_t> corrections_compress_fastpfor;

        sdsl::int_vector<64> corrections_none; // corrections for saving, each value <= (epsilon / 4)
        sdsl::int_vector<64> corrections_exception; // corrections for saving, each (epsilon / 4) < value <= epsilon
        std::vector<Correction_Value> corrections_vector; // corrections for decode

        uint64_t segments_size;
        uint8_t bpc_first, bpc_covered, bpc_intercept, bpc_slope_significand, bpc_slope_exponent;

        /// Sentinel value to avoid bounds checking.
        static constexpr K sentinel = std::numeric_limits<K>::has_infinity ? std::numeric_limits<K>::infinity() : std::numeric_limits<K>::max();

        using position_type = typename std::conditional_t<sizeof(K) <= 4, uint32_t, uint64_t>;

        using canonical_segment = typename internal::OptimalPiecewiseLinearModel<position_type, K>::CanonicalSegment;

        uint64_t errorPointCount;

        PGM(): n(0), errorPointCount(0) {}

        explicit PGM(const std::vector<K>& data) : PGM(data.begin(), data.end()) {}

        template <typename RandomIt>
        PGM(RandomIt begin, RandomIt end):
            n(std::distance(begin, end)),
            first_pos(n ? *begin : K(0)),
            levels_offsets(),
            errorPointCount(0){
            Epsilon_Data = Epsilon;
            corrections_vector.resize(n);

            if (residual_compress_type == "none"){
                bit_compress = 0;
                epsilon_compress = Epsilon;
                bpc_compress = BIT_WIDTH(Epsilon);
            } else {
                bit_compress = 1;
            }

            if (bit_compress == 0) {
                corrections_none.resize(((n * bpc_compress + 63) >> 6)); // set 3 to avoid insEufficient memory
                signs_none.resize(n);
            }

            build(begin, end, Epsilon,  levels_offsets, errorPointCount, corrections_none, corrections_exception,  signs_none, signs_exception);

            if (residual_compress_type == "fastpfor")
                corrections_compress_fastpfor = compress_residuals_fastpfor(corrections_vector);

            segments_compress();

        }

        template <typename RandomIt>
        void build(RandomIt begin, RandomIt end,
                          uint64_t epsilon,
                          std::vector<uint64_t>& levels_offsets,
                          uint64_t& errorPointCount,
                          sdsl::int_vector<64>& corrections_compress,
                          sdsl::int_vector<64>& corrections_exception,
                          sdsl::bit_vector& signs_compress,
                          sdsl::bit_vector& signs_exception) {

            auto n = (uint64_t) std::distance(begin, end);
            if (n == 0) return;

            levels_offsets.push_back(0);

            if (*std::prev(--end) == sentinel)
                throw std::invalid_argument("The value " + std::to_string(sentinel) + " is reserved as a sentinel.");

            std::vector<canonical_segment> canonical_segments;
            canonical_segments.reserve(epsilon > 0 ? n / (epsilon * epsilon) : n / 8);

            auto in_fun = [begin](auto i) { return std::pair<position_type, K>(i, begin[i]); };
            auto out_fun = [&canonical_segments](auto cs) { canonical_segments.push_back(cs); };
            auto n_segments = internal::make_segmentation_par(n, epsilon, in_fun, out_fun);
            auto last_n = n_segments;

            for (auto it = canonical_segments.begin(); it < canonical_segments.end(); ++it) {
                auto i = it->get_first_x();
                auto j = std::next(it) != canonical_segments.end() ? std::next(it)->get_first_x() : n;
                build_segments(*it, n, i, j, begin, errorPointCount, corrections_compress.data(), corrections_exception.data(), signs_compress, signs_exception); // build the segment
            }
            levels_offsets.push_back(seg_first.size());
        }

        template <typename RandomIt>
        void build_segments(const canonical_segment& cs,
            uint64_t n, uint64_t i, uint64_t j,
            RandomIt data, uint64_t& errorPointCount,
            uint64_t* corrections_compress_ptr, uint64_t* corrections_exception_ptr,
            sdsl::bit_vector& signs_compress_ptr, sdsl::bit_vector& signs_exception_ptr) {

            uint32_t first = cs.get_first_x();
            if (first == n) return;

            seg_first.push_back(first);

            auto [cs_significand, cs_exponent, cs_intercept] = cs.get_fixed_point_segment(first, j - i + 1); // fixed point slope and intercept

            if (first == n - 1) {
                seg_intercept.push_back(data[first]);
                seg_slope_exponent.push_back(0);
                seg_slope_significand.push_back(0);
                seg_covered.push_back(1);

                if (bit_compress > 0) {
                    corrections_vector[first] = 0;
                } else {
                    set_correction(corrections_compress_ptr, first, 0, 0, signs_compress_ptr);
                }
                return;
            }

            seg_slope_exponent.push_back(cs_exponent);
            seg_slope_significand.push_back(cs_significand);
            seg_intercept.push_back(cs_intercept);
            seg_covered.push_back(j - i);


            if (residual_compress_type == "fastpfor") {
                int64_t last_correction = 0;
                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    int64_t error_diff = error - last_correction;
                    corrections_vector[p] = error_diff;
                    last_correction = error;
                }
            } else if (residual_compress_type == "none") {
                for (Covered_Value p = first; p < j; p++) {
                    int64_t error = static_cast<int64_t> (data[p]) - seg_approximate(p, first, cs_exponent, cs_significand, cs_intercept);
                    uint8_t sign_value = error >= 0 ? 0 : 1;
                    error = std::abs(error);
                    if (error <= Epsilon)
                        set_correction(corrections_compress_ptr, p, error, sign_value, signs_compress_ptr);
                    else if (error == Epsilon + 1 && sign_value == 0) // use -0 to represent error = epsilon
                        set_correction(corrections_compress_ptr, p, 0, 1, signs_compress_ptr);
                    else
                        std::cerr << "Error: error = -epsilon: " << error << " " << sign_value << std::endl;
                }
            } else {
                throw std::runtime_error("residual_compress_type not recognised");
            }
        }

        void segments_compress() {
            uint32_t max_fist = seg_first.back();
            bpc_first = BIT_WIDTH(max_fist);
            auto max_iter_1 = std::max_element(seg_covered.begin(), seg_covered.end());
            bpc_covered = BIT_WIDTH(*max_iter_1);
            auto max_iter_2 = std::max_element(seg_intercept.begin(), seg_intercept.end());
            bpc_intercept = BIT_WIDTH(*max_iter_2);
            auto max_iter_3 = std::max_element(seg_slope_significand.begin(), seg_slope_significand.end());
            bpc_slope_significand = BIT_WIDTH(*max_iter_3);
            auto max_iter_4 = std::max_element(seg_slope_exponent.begin(), seg_slope_exponent.end());
            bpc_slope_exponent = BIT_WIDTH(*max_iter_4);

            segments_size = seg_first.size();
            seg_first_compress.resize((segments_size * bpc_first + 63) >> 6);
            seg_covered_compress.resize((segments_size * bpc_covered + 63) >> 6);
            seg_intercept_compress.resize((segments_size * bpc_intercept + 63) >> 6);
            seg_slope_significand_compress.resize((segments_size * bpc_slope_significand + 63) >> 6);
            seg_slope_exponent_compress.resize((segments_size * bpc_slope_exponent + 63) >> 6);

            for (int i = 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
                set_segment_compress(seg_first_compress.data(), i, first, bpc_first);
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
            seg_first.clear();
            seg_covered.clear();
            seg_intercept.clear();
            seg_slope_significand.clear();
            seg_slope_exponent.clear();
            corrections_vector.clear();
            std::vector<K>().swap(seg_first);
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Intercept_Value>().swap(seg_intercept);
            std::vector<int64_t>().swap(seg_slope_significand);
            std::vector<uint8_t>().swap(seg_slope_exponent);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        int64_t seg_approximate(uint32_t i, uint32_t first, uint8_t exponent, int64_t significand, int32_t intercept) const {
            return (int64_t(significand * (i - first)) >> exponent) + intercept;
        }

        uint64_t segments_count() const { return levels_offsets[1] - 1; }

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

        uint64_t segment_size_in_bytes() const {
            // return ((seg_covered_compress.bit_size() + seg_intercept_compress.bit_size() + seg_slope_significand_compress.bit_size() + seg_slope_exponent_compress.bit_size()) / CHAR_BIT) + sizeof(Epsilon_Data);
            return ((seg_first_compress.bit_size() + seg_intercept_compress.bit_size() + seg_slope_significand_compress.bit_size() + seg_slope_exponent_compress.bit_size()) / CHAR_BIT) + sizeof(Epsilon_Data);
        }

        uint64_t corrections_size_in_bytes() const {
            if (bit_compress == 0)
                return (corrections_none.bit_size()) / CHAR_BIT;
            return corrections_compress_fastpfor.size() * 4;
        }

        uint64_t signs_size_in_bytes() const {
            if (bit_compress == 0)
                return (signs_none.bit_size()) / CHAR_BIT;
            return 0;
        }

        void report_residual_random_segment(std::ofstream &file, int seg_idx) {
            if (bit_compress > 0) {
                corrections_exception_num_read = 0;
                for (int i = 0; i < segments_size; i++) {
                    if (i == seg_idx) {
                        auto first = seg_first[i];
                        auto covered = seg_covered[i];
                        int32_t last_correction = 0;
                        file << covered << std::endl;
                        for (int j = first; j < covered + first; j++) {
                            auto error1 = get_correction_compress(corrections_none.data(), corrections_exception.data(), j, signs_none, signs_exception);
                            auto error2 = last_correction + error1;
                            last_correction = error2;
                            auto flag = std::abs(error1) <= epsilon_compress ? 0 : 1;
                            file << error2 << "\t" << error1 << "\t"<< flag << std::endl;
                        }
                    }
                }
            } else {
                throw std::invalid_argument("Not support uncompress residuals");
            }
        }

        void report_residual_list(std::ofstream &file) {
            if (bit_compress > 0) {
                corrections_exception_num_read = 0;
                file << n << std::endl;
                for (int i = 0; i < segments_size; i++) {
                    auto first = seg_first[i];
                    auto covered = seg_covered[i];
                    int32_t last_correction = 0;
                    // file << covered << std::endl;
                    for (int j = first; j < covered + first; j++) {
                        // auto error1 = get_correction(corrections_compress.data(), corrections_exception.data(), j, signs_compress, signs_exception);
                        auto delta_residual = corrections_vector[j];
                        auto true_residual = last_correction + delta_residual;
                        last_correction = true_residual;
                        // auto flag = std::abs(error1) <= epsilon_compress ? 0 : 1;
                        file << true_residual << " " << delta_residual << " " << zigzag(delta_residual) << std::endl;
                    }
                }
            } else {
                throw std::invalid_argument("Not support uncompress residuals");
            }
        }

        void segment_init() {
            seg_first.resize(segments_size);
            seg_covered.resize(segments_size);
            seg_intercept.resize(segments_size);
            seg_slope_significand.resize(segments_size);
            seg_slope_exponent.resize(segments_size);

            for (int i = 0; i < segments_size; i++) {
                seg_first[i] = get_segment_compress(seg_first_compress.data(), i, bpc_first);
                seg_covered[i] = get_segment_compress(seg_covered_compress.data(), i, bpc_covered);
                seg_intercept[i] = get_segment_compress(seg_intercept_compress.data(), i, bpc_intercept);
                seg_slope_significand[i] = get_segment_compress(seg_slope_significand_compress.data(), i, bpc_slope_significand);
                seg_slope_exponent[i] = get_segment_compress(seg_slope_exponent_compress.data(), i, bpc_slope_exponent);
            }
        }

        void residual_init() {
            if (residual_compress_type == "fastpfor") {
                std::vector<uint32_t> uncompressed_output(n);
                corrections_vector.resize(n);
                decompress_residuals_fastpfor(corrections_compress_fastpfor, n, uncompressed_output, corrections_vector);
            } else if (residual_compress_type == "none") { // for uncompressed corrections, we need to change them to delta type for decode
                corrections_vector.resize(n);
                for(int i = 0; i < seg_first.size(); i++) {
                    auto covered = seg_covered[i];
                    auto first = seg_first[i];
                    int32_t last_correction = 0;
                    for (Covered_Value j = first; j < first + covered; ++j) {
                        int32_t correction_varify = get_correction_uncompress(corrections_none.data(), j, signs_none);
                        corrections_vector[j] = correction_varify - last_correction;
                        last_correction = correction_varify;
                    }
                }
            } else {
                throw std::runtime_error("Not support uncompress residuals");
            }
        }

        void normal_init(){
            segment_init();
            residual_init();
        }

        void normal_clean() {
            seg_first.clear();
            seg_covered.clear();
            seg_intercept.clear();
            seg_slope_significand.clear();
            seg_slope_exponent.clear();
            corrections_vector.clear();
            std::vector<K>().swap(seg_first);
            std::vector<Covered_Value>().swap(seg_covered);
            std::vector<Intercept_Value>().swap(seg_intercept);
            std::vector<int64_t>().swap(seg_slope_significand);
            std::vector<uint8_t>().swap(seg_slope_exponent);
            std::vector<Correction_Value>().swap(corrections_vector);
        }

        std::vector<K> normal_decode_rice() {
            std::vector<K> output;
            output.resize(n);
            K* output_pointer = output.data();

            uint32_t max_covered = 0;
            for (int i= 0; i < seg_first.size(); i++) {
                max_covered = seg_covered[i] > max_covered ? seg_covered[i] : max_covered;
            }

            const uint32_t residuals_block_size = ((max_covered + block_size_rice - 1) / block_size_rice) * block_size_rice;
            std::vector<int> residuals_block_vector(residuals_block_size);
            int *residuals_block = residuals_block_vector.data();

            uint32_t block_remain = 0;
            uint32_t block_pointer = 0;

            RiceDecoder rd(corrections_compress_rice.data(), corrections_compress_rice.size());

            for (int i= 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
                auto covered = seg_covered[i];
                auto significand = seg_slope_significand[i];
                auto exponent = seg_slope_exponent[i];
                auto intercept = seg_intercept[i];


                int last_correction = intercept;

                int j = 0;
                while (j < covered) {
                    if (block_remain == 0) {
                        block_pointer = 0;
                        uint32_t decode_length = rd.decompress_residuals_block(residuals_block);
                        block_remain = decode_length;

                        while (block_remain < residuals_block_size && decode_length > 0) {
                            decode_length = rd.decompress_residuals_block(residuals_block + block_remain);
                            block_remain += decode_length;
                        }

                        if (block_remain == 0) {
                            throw std::runtime_error("Unexpected end of residual stream");
                        }
                    }
                    // 本次最多消费这么多：块剩余 vs. 段剩余
                    int take = std::min<int>(block_remain, covered - j);

                    for (int pos = j; pos < j + take; ++pos) {
                        last_correction += residuals_block[block_pointer++];
                        *output_pointer++ = ((pos * significand) >> exponent) + last_correction;
                    }

                    j += take;
                    block_remain -= take;
                }
            }
            return output;
        }

        std::vector<K> normal_decode_fastpfor() {
            std::vector<K> output;
            output.resize(n);
            uint32_t pointer = 0;

            std::vector<uint32_t> uncompressed_output;
            uncompressed_output.resize(n);

            corrections_vector.resize(n);

            decompress_residuals_fastpfor(corrections_compress_fastpfor, n, uncompressed_output, corrections_vector);

            for (int i= 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
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

        std::vector<K> normal_decode() {
            std::vector<K> output;
            output.resize(n);
            uint32_t pointer = 0;

            for (int i= 0; i < seg_first.size(); i++) {
                auto first = seg_first[i];
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

        void set_correction(uint64_t* corrections_compress_val, uint64_t i, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_compress_ptr) {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            sdsl::bits::write_int(corrections_compress_val + (idx >> 6u), value, idx & 0x3f, bpc_compress);
            signs_compress_ptr[i] = sign_value & 0x01;
        }

        void set_correction_exception(uint64_t* corrections_exception_val, uint64_t value, uint8_t sign_value, sdsl::bit_vector& signs_exception_ptr) {
            if (value > Epsilon_Data && Epsilon_Data == 1) {
                std::cerr << "Correction Exception Value: " << value << std::endl;
            }
            uint64_t idx = get_correction_bit_offset(corrections_exception_num_write, bpc_exception);
            sdsl::bits::write_int(corrections_exception_val + (idx >> 6u), value, idx & 0x3f, bpc_exception);
            signs_exception_ptr[corrections_exception_num_write++] = sign_value & 0x01;
        }

        void set_correction_compress(uint64_t* corrections_compress_val, uint64_t* corrections_exception_val, uint64_t i, int64_t value, sdsl::bit_vector& signs_compress_ptr, sdsl::bit_vector& signs_exception_ptr) {
            uint64_t abs_value = std::abs(value);
            if (abs_value <= epsilon_compress) {
                set_correction(corrections_compress_val, i, abs_value, value >= 0 ? 0 : 1, signs_compress_ptr);
            } else { // if is unexpected value
                set_correction_exception(corrections_exception_val, abs_value, value >= 0 ? 0 : 1, signs_exception_ptr); // store the exception value
                set_correction(corrections_compress_val, i, 0, 1, signs_compress_ptr); // special sign -0 to record it's an exception value
            }
        }

        int64_t get_correction_compress(const uint64_t* corrections_compress_val, const uint64_t* corrections_exception_val, int64_t i, const sdsl::bit_vector& signs_compress_ptr, const sdsl::bit_vector& signs_exception_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_compress);
            if (correction == 0 && signs_compress_ptr[i] == 1) {
                idx = get_correction_bit_offset(corrections_exception_num_read, bpc_exception);
                correction = sdsl::bits::read_int(corrections_exception_val + (idx >> 6u), idx & 0x3f, bpc_exception);
                return signs_exception_ptr[corrections_exception_num_read++] == 0 ? correction : -correction;
            }
            return signs_compress_ptr[i] == 0 ? correction : -correction;
        }

        int64_t get_correction_uncompress(const uint64_t* corrections_compress_val, int64_t i, const sdsl::bit_vector& signs_compress_ptr)  {
            uint64_t idx = get_correction_bit_offset(i, bpc_compress);
            uint64_t correction = sdsl::bits::read_int(corrections_compress_val + (idx >> 6u), idx & 0x3f, bpc_compress);
            if (correction == 0 && signs_compress_ptr[i] == 1)
                return Epsilon + 1;
            return signs_compress_ptr[i] == 0 ? correction : -correction;
        }

     };

}
