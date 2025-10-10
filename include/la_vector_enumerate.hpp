#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <config.hpp>
#include  <variant>
#include <lico_index.hpp>
#include <tools.hpp>


namespace la_vector_sequence {
    using LICOVariant = std::variant<
        lico::LICO<uint32_t, 3>,
        lico::LICO<uint32_t, 7>,
        lico::LICO<uint32_t, 15>,
        lico::LICO<uint32_t, 31>,
        lico::LICO<uint32_t, 63>,
        lico::LICO<uint32_t, 127>,
        lico::LICO<uint32_t, 255>,
        lico::LICO<uint32_t, 511>,
        lico::LICO<uint32_t, 1023>,
        lico::LICO<uint32_t, 2047>,
        lico::LICO<uint32_t, 4095>,
        lico::LICO<uint32_t, 8191>,
        lico::LICO<uint32_t, 16383>,
        lico::LICO<uint32_t, 32767>,
        lico::LICO<uint32_t, 65535>,
        lico::LICO<uint32_t, 131071>,
        lico::LICO<uint32_t, 262143>,
        lico::LICO<uint32_t, 524287>,
        lico::LICO<uint32_t, 1048575>>;


    template <typename K>
    class la_vector_enumerator{

        typedef int64_t Simd_Value;
        typedef uint64_t Slope_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Intercept_Value;
        typedef uint32_t Covered_Value;

        public:

        struct segment{
            Covered_Value first;
            Intercept_Value intercept; // 32 bits
            uint8_t slope_exponent;
            Slope_Value slope_significand;
            Covered_Value covered; // covered
            inline segment(Covered_Value first, Intercept_Value intercept, uint8_t slope_exponent, Slope_Value slope_significand, Covered_Value covered) :
                    first(first), intercept(intercept), slope_exponent(slope_exponent), slope_significand(slope_significand), covered(covered) {}
        };

        uint64_t n;

        std::vector<segment> segments;

        std::vector<K> segment_first_values; // size == segments.size()

        std::vector<K> segment_max;

        std::vector<uint64_t> parted_sizes;

        std::vector<Correction_Value> corrections_vector; // corrections for decode


        std::vector<uint32_t> block_sizes;

        void load_block_size(std::vector<uint32_t> all_block_sizes) {
            this -> block_sizes = std::move(all_block_sizes);
        }

        void load_residuals(uint64_t data_size, std::vector<Correction_Value> corrections_vector) {
            this -> n = data_size;
            this -> corrections_vector = std::move(corrections_vector);
        }

        // used for Query Test
        K current_value = INT_MAX;
        K next_first_value = INT_MAX;
        Correction_Value current_correction = 0;
        Covered_Value current_pos = 0;
        Covered_Value current_segment = 0;
        uint32_t total_segment_size = 0;
        std::vector<K> current_value_vector;

        void warm_up() {
            for (auto k = 0; k < 5; k++) {
                 for (auto i = 0;i < n; i++)
                     current_value_vector[i] = 0;
            }
        }

        void query_init(const std::string decode_type, const std::string query_type) {
            if (decode_type == "normal") {
                if (query_type == "intersection") {
                    current_pos = 0;
                    current_segment = 0;
                    total_skip = 0;
                    // current_value_vector.resize(n);
                } else if (query_type == "union") {
                    current_pos = 0;
                    current_segment = 0;
                    current_value = INT_MAX - 1;
                    total_skip = 0;
                    // current_value_vector.resize(n);
                }
            } else if (decode_type == "simd") {
                current_pos = 0;
                current_segment = 0;
                // current_value_vector.resize(n);
                // simd_init();
            }
        }

        void decode_query(K* output, const std::string decode_type) {
            normal_decode(output);
        }

        long double total_skip = 0;

        void build_segment_first_values() {
            segment_first_values.clear();
            segment_first_values.resize(segments.size() + 1);
            K posi = 0;
            auto start = std::chrono::high_resolution_clock::now();
            for (auto &s : segments) {
                segment_first_values[posi++] = s.intercept + corrections_vector[s.first];
            }
            segment_first_values[posi] = INT_MAX;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        uint64_t total_duration = 0;

        K* next_pointer;

        void next_init() {
            current_pos = 0;
            current_value_vector.resize(n + 1);

            for (Covered_Value j = 0; j < n; j++) // warm up
                current_value_vector[j] = j + 1;

            normal_decode(current_value_vector.data());
            current_value_vector[n] = INT_MAX;
            next_pointer = current_value_vector.data();
            current_value = *next_pointer;
        }

        K docid() {
            return current_value;
        }

        void next() {
            current_value = *++next_pointer;
        }

        void next_geq_init() {
            residuals_decode();
            build_segment_first_values();

            total_segment_size = segments.size();

            current_segment = 0;
            if (segments.size() > 1)
                this -> next_first_value = segment_first_values[1];
            else
                this -> next_first_value = INT_MAX;

            this -> current_correction = segments[0].intercept;
            this -> current_value = segment_first_values[0];

        }

        K nextgeq(K posting_value) {
            if (current_value >= posting_value) return current_value;

            while (posting_value >= next_first_value) {
                ++current_segment;
                if (current_segment >= total_segment_size) { current_pos = n; return INT_MAX; }
                auto &cur_seg = segments[current_segment];
                current_pos = cur_seg.first;
                current_correction = cur_seg.intercept;
                next_first_value = segment_first_values[current_segment + 1];
            }

            while (true){
                auto &seg = segments[current_segment];

                Covered_Value j_base = current_pos > seg.first ? (current_pos - seg.first) : 0;
                Slope_Value accum = (static_cast<Slope_Value>(j_base) * seg.slope_significand);
                const Correction_Value* corrections_pointer = corrections_vector.data() + j_base + seg.first;
                for (Covered_Value j = j_base; j < seg.covered; ++j) {
                    current_correction += *corrections_pointer++;
                    current_value = static_cast<Correction_Value>(accum >> seg.slope_exponent) + current_correction;
                    if (current_value >= posting_value) {
                        current_pos = seg.first + j + 1;
                        return current_value;
                    }
                    accum += seg.slope_significand;
                }

                ++current_segment;
                if (current_segment >= total_segment_size) { current_pos = n; return INT_MAX; }
                auto &cur_seg = segments[current_segment];
                current_pos = cur_seg.first;
                current_correction = cur_seg.intercept;
                next_first_value = segment_first_values[current_segment + 1];
            }
        }



        void normal_decode(K* output) {
            Correction_Value* correction_pointer = corrections_vector.data();
            auto start = std::chrono::high_resolution_clock::now();
            const auto end_iter = segments.end();
            for (auto it = segments.begin(); it != end_iter; ++it) {
                const auto& seg = *it;
                Slope_Value significand =0;
                Correction_Value last_correction = seg.intercept;
                for (Covered_Value j = 0; j < seg.covered; ++j) {
                    last_correction = last_correction + *correction_pointer++;
                    *output++ = static_cast<Correction_Value> (significand >> seg.slope_exponent) + last_correction;
                    significand += seg.slope_significand;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        void free_memory(std::string decode_type = "simd") {
            std::vector<K> ().swap(current_value_vector);
            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);
        }

    
    };

    template <typename K>
    static la_vector_enumerator<K> create_enumerator_from_single_index(LICOVariant &index) {
        la_vector_enumerator<K> enumerator;

        uint32_t last_first = 0;

        std::visit([&enumerator, &last_first](auto &idx) {
            idx.segment_init();

            for (int i = 0; i < idx.segments_size; ++i) {
                auto first = last_first;
                auto intercept = idx.seg_intercept[i];
                auto slope_exponent = idx.seg_slope_exponent[i];
                auto slope_significand = idx.seg_slope_significand[i];
                auto covered = idx.seg_covered[i];
                last_first += covered;
                enumerator.segments.emplace_back(first, intercept, slope_exponent, slope_significand, covered);
            }

                idx.residual_init();
                enumerator.load_residuals(idx.n, idx.corrections_vector);
            enumerator.load_block_size(std::vector<uint32_t>{idx.segments_size});
            idx.free_memory();
        }, index);

        return enumerator;
    }

    template <typename K, size_t epsilon = 64> // K is uint32_t or uint64_t
    class la_vector_decoder{

    public:
        uint64_t data_size = 0;
        uint64_t total_list = 0;

        uint64_t total_decode_time = 0;
        uint64_t max_decode_time = 0;
        uint64_t min_decode_time = UINT64_MAX - 1;

        uint64_t total_calculated = 0;
        uint64_t total_calculated_add = 0;
        uint64_t total_conversion_time = 0;
        uint64_t total_unequal = 0;
        std::vector<std::string> perf_header;
        std::vector<std::vector<double>> perf_data;

        void decode_test(la_vector_enumerator <K> &index, std::string &decode_type, bool warm_up = true) {
            total_list++;

            std::vector<K> result2(index.n);
            if (warm_up) {
                for (int warm_time = 0; warm_time < 5; warm_time++) {
                    for (int i = 0; i < index.n; i ++) {
                        result2[i] = 0;
                    }
                }
            }
            index.normal_decode(result2.data());

            std::vector<K> ().swap(result2);

            if (index.total_duration > max_decode_time)
                max_decode_time = index.total_duration;
            if (index.total_duration < min_decode_time)
                min_decode_time = index.total_duration;

            index.free_memory(decode_type);
            data_size += index.n;
            total_decode_time += index.total_duration;
        }

        void result_statistic(std::string &decode_type) {
            std::cerr << "Total list: " << total_list << std::endl;
            if (decode_type == "simd" || decode_type == "simd_simple") {
                std::cerr << "Decode time 1 simd, average: " << total_decode_time / total_list  << ", max: " << max_decode_time << ", min: " << min_decode_time << " microseconds" <<std::endl;
                std::cerr << "Total calculated: " << total_calculated << " (" << double(total_calculated) / data_size << ") , Total calculated add: " << total_calculated_add << ", Total Conversion time: " << total_conversion_time << ", Total Unequal: " << total_unequal << std::endl;
            }
            if (decode_type == "normal") {
                std::cerr << "Decode time 2 normal, average: " << total_decode_time / total_list  << ", max: " << max_decode_time << ", min: " << min_decode_time << " microseconds" <<std::endl;
            }

            std::cerr << "Decode per integer: " << static_cast<long double> (total_decode_time) / data_size * 1000.0 << " nanoseconds" << std::endl;

        }

        void test_model(const std::string input_basename, std::string decode_type) {
            if (input_basename.empty() || input_basename.back() != '/') {
                std::cerr << "Error: input_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Load index from: " << input_basename << std::endl;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                return;
            }


            uint64_t data_size_tmp = 0;
            in_header.read(reinterpret_cast<char*>(&data_size_tmp), sizeof(data_size_tmp));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            for (K i = 0; i < index_num; ++i) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);

                uint32_t Epsilon_Data;

                while (in.peek() != EOF) {
                    LICOVariant variant_index;
                    in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));
                    switch (Epsilon_Data) {
                        case 3: variant_index = lico::LICO<K, 3>(); break;
                        case 7: variant_index = lico::LICO<K, 7>(); break;
                        case 15: variant_index = lico::LICO<K, 15>(); break;
                        case 31: variant_index = lico::LICO<K, 31>(); break;
                        case 63: variant_index = lico::LICO<K, 63>(); break;
                        case 127: variant_index = lico::LICO<K, 127>(); break;
                        case 255: variant_index = lico::LICO<K, 255>(); break;
                        case 511: variant_index = lico::LICO<K, 511>(); break;
                        case 1023: variant_index = lico::LICO<K, 1023>(); break;
                        case 2047: variant_index = lico::LICO<K, 2047>(); break;
                        case 4095: variant_index = lico::LICO<K, 4095>(); break;
                        case 8191: variant_index = lico::LICO<K, 8191>(); break;
                        case 16383: variant_index = lico::LICO<K, 16383>(); break;
                        case 32767: variant_index = lico::LICO<K, 32767>(); break;
                        case 65535: variant_index = lico::LICO<K, 65535>(); break;
                        case 131071: variant_index = lico::LICO<K, 131071>(); break;
                        case 262143: variant_index = lico::LICO<K, 262143>(); break;
                        case 524287: variant_index = lico::LICO<K, 524287>(); break;
                        case 1048575: variant_index = lico::LICO<K, 1048575>(); break;
                        default: std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl; break;
                    }

                    std::visit([&in, Epsilon_Data](auto& index) {
                        index.Epsilon_Data = Epsilon_Data;
                        read_index_data(in, index);
                    }, variant_index);

                    la_vector_enumerator <K> enumerator_tmp = create_enumerator_from_single_index <K> (variant_index);

                    decode_test(enumerator_tmp, decode_type);
                }

                in.close();
            }


            in_header.close();
            result_statistic(decode_type);
        }

    };

    template <typename K, uint64_t epsilon = 64> // K is uint32_t or uint64_t, Floating is unused
    class lico_querier{
    private:
        uint64_t data_size = 0;
        uint32_t query_num = 0;
        uint64_t query_time = 0;
        std::string input_basename = "";

        std::vector<std::vector<uint32_t>> read_query(const std::string& filename) {
            std::vector<std::vector<uint32_t>> idLists;
            std::ifstream file(filename);
            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << filename << std::endl;
                return idLists; // Return an empty vector if the file could not be opened.
            }
            std::string line;
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::vector<uint32_t> ids;
                uint32_t id;
                while (iss >> id) {                 // Extract uint32_t from the line until no more can be found.
                    ids.push_back(id);
                }
                remove_duplicate_terms(ids);
                idLists.push_back(ids);
            }
            query_num = idLists[0].size();
            std::cout << "Total query sequences: " << idLists.size() << " Query num: " << (query_num > 5 ? 5 : query_num) << std::endl;
            file.close();
            return idLists;
        }

        std::vector<la_vector_enumerator <K>> load_model(std::vector<uint32_t> idx_list) {
            if (input_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                throw std::runtime_error("file format error");
            }

            std::vector<la_vector_enumerator <K>> index_sequences;

            std::ifstream in_header(input_basename + "idx.size", std::ios::binary);
            if (!in_header) {
                std::cerr << "Error: Cannot open idx.size for reading." << std::endl;
                throw std::runtime_error("File open error");
            }

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            for (auto i : idx_list) {
                if (i >= index_num) {
                    std::cerr << "Index count out of range: " << i << "Total index: " << index_num << std::endl;
                    return {};
                }
                std::string filename = input_basename + std::to_string(i) + ".idx";

                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    throw std::runtime_error("File open error");
                }

                uint32_t Epsilon_Data;
                in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                LICOVariant variant_index;
                switch (Epsilon_Data) {
                    case 3: variant_index = lico::LICO<K, 3>(); break;
                    case 7: variant_index = lico::LICO<K, 7>(); break;
                    case 15: variant_index = lico::LICO<K, 15>(); break;
                    case 31: variant_index = lico::LICO<K, 31>(); break;
                    case 63: variant_index = lico::LICO<K, 63>(); break;
                    case 127: variant_index = lico::LICO<K, 127>(); break;
                    case 255: variant_index = lico::LICO<K, 255>(); break;
                    case 511: variant_index = lico::LICO<K, 511>(); break;
                    case 1023: variant_index = lico::LICO<K, 1023>(); break;
                    case 2047: variant_index = lico::LICO<K, 2047>(); break;
                    case 4095: variant_index = lico::LICO<K, 4095>(); break;
                    case 8191: variant_index = lico::LICO<K, 8191>(); break;
                    case 16383: variant_index = lico::LICO<K, 16383>(); break;
                    default:
                        std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl;
                        continue;
                }

                std::visit([&in, Epsilon_Data](auto& index) {
                    index.Epsilon_Data = Epsilon_Data;
                    read_index_data(in, index);
                }, variant_index);

                in.close();

                la_vector_enumerator <K> enumerator_tmp = create_enumerator_from_single_index<K>(variant_index);

                index_sequences.emplace_back(enumerator_tmp);
            }

            in_header.close();

            return index_sequences;
        }

        static uint32_t intersect_u32_normal(uint32_t const *a, uint32_t const *b,uint32_t a_length, uint32_t b_length, uint32_t *out) {
            uint32_t* const intersect_start = out;

            const uint32_t* a_end = a + a_length;
            const uint32_t* b_end = b + b_length;

            while (a < a_end && b < b_end) {
                if (*a < *b) {
                    ++a;
                } else if (*a > *b) {
                    ++b;
                } else {
                    // Found a match
                    *out++ = *a;
                    ++a;
                    ++b;
                }
            }
            return out - intersect_start;
        }

        static uint32_t intersect_u32_simd_basic(uint32_t const* shorter, uint32_t const* longer, uint32_t shorter_length, uint32_t longer_length, uint32_t * out) {
            uint32_t intersection_count = 0;
            uint32_t shorter_idx = 0, longer_idx = 0;
            uint32_t longer_load_size;
            __mmask16 longer_mask;

            while (shorter_idx < shorter_length && longer_idx < longer_length) {
                // Load `shorter_member` and broadcast it to shorter vector, load `longer_members_vec` from memory.
                uint32_t longer_remaining = longer_length - longer_idx;
                uint32_t shorter_member = shorter[shorter_idx];
                __m512i shorter_member_vec = _mm512_set1_epi32(*(int*)&shorter_member);
                __m512i longer_members_vec;
                if (longer_remaining < 16) {
                    longer_load_size = longer_remaining;
                    longer_mask = (__mmask16)_bzhi_u32(0xFFFF, longer_remaining);
                } else {
                    longer_load_size = 16;
                    longer_mask = 0xFFFF;
                }
                longer_members_vec = _mm512_maskz_loadu_epi32(longer_mask, (__m512i const*)(longer + longer_idx));

                // Compare `shorter_member` with each element in `longer_members_vec`,
                // and jump to the position of the match. There can be only one match at most!
                __mmask16 equal_mask = _mm512_mask_cmpeq_epu32_mask(longer_mask, shorter_member_vec, longer_members_vec);
                bool equal_count = equal_mask != 0;
                if (equal_count) {
                    out[intersection_count++] = shorter_member;
                }

                // When comparing a scalar against a sorted array, we can find three types of elements:
                // - entries that scalar is greater than,
                // - entries that scalar is equal to,
                // - entries that scalar is less than,
                // ... in that order! Any of them can be an empty set.
                __mmask16 greater_mask = _mm512_mask_cmplt_epu32_mask(longer_mask, longer_members_vec, shorter_member_vec);
                uint32_t greater_count = _mm_popcnt_u32(greater_mask);
                uint32_t smaller_exists = longer_load_size > greater_count - equal_count;

                // Advance the first array:
                // - to the next element, if a match was found,
                // - to the next element, if the current element is smaller than any elements in the second array.
                shorter_idx += equal_count | smaller_exists;
                // Advance the second array:
                // - to the next element after match, if a match was found,
                // - to the first element that is greater than the current element in the first array, if no match was found.
                longer_idx += greater_count + equal_count;

                // At any given cycle, take one entry from shorter array and compare it with multiple from the longer array.
                // For that, we need to swap the arrays if necessary.
                if ((shorter_length - shorter_idx) > (longer_length - longer_idx)) {
                    uint32_t const* temp_array = shorter;
                    shorter = longer, longer = temp_array;
                    uint32_t temp_length = shorter_length;
                    shorter_length = longer_length, longer_length = temp_length;
                    uint32_t temp_idx = shorter_idx;
                    shorter_idx = longer_idx, longer_idx = temp_idx;
                }
            }
            return intersection_count;
        }

        static uint32_t intersect_u32_simd(uint32_t const* a, uint32_t const* b, uint32_t a_length, uint32_t b_length, uint32_t *out) {
            uint32_t const* const a_end = a + a_length;
            uint32_t const* const b_end = b + b_length;
            uint32_t c = 0;
            union vec_t {
                __m512i zmm;
                uint32_t u32[16];
            } a_vec, b_vec;

            while (a + 16 < a_end && b + 16 < b_end) {
                a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
                b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);

                // Intersecting registers with `_mm512_2intersect_epi16_mask` involves a lot of shuffling
                // and comparisons, so we want to avoid it if the slices don't overlap at all
                uint32_t a_min;
                uint32_t a_max = a_vec.u32[15];
                uint32_t b_min = b_vec.u32[0];
                uint32_t b_max = b_vec.u32[15];

                // If the slices don't overlap, advance the appropriate pointer
                while (a_max < b_min && a + 32 < a_end) {
                    a += 16;
                    a_vec.zmm = _mm512_loadu_si512((__m512i const*)a);
                    a_max = a_vec.u32[15];
                }
                a_min = a_vec.u32[0];
                while (b_max < a_min && b + 32 < b_end) {
                    b += 16;
                    b_vec.zmm = _mm512_loadu_si512((__m512i const*)b);
                    b_max = b_vec.u32[15];
                }
                b_min = b_vec.u32[0];

                // Now we are likely to have some overlap, so we can intersect the registers
                __mmask16 a_matches = _mm512_2intersect_epi32_mask(a_vec.zmm, b_vec.zmm);
                _mm512_mask_compressstoreu_epi32(out + c, a_matches, a_vec.zmm);


                c += _mm_popcnt_u32(a_matches); // The `_popcnt32` symbol isn't recognized by MSVC

                // Determine the number of entries to skip in each array, by comparing
                // every element in the vector with the last (largest) element in the other array
                __m512i a_last_broadcasted = _mm512_set1_epi32(*(int const*)&a_max);
                __m512i b_last_broadcasted = _mm512_set1_epi32(*(int const*)&b_max);
                __mmask16 a_step_mask = _mm512_cmple_epu32_mask(a_vec.zmm, b_last_broadcasted);
                __mmask16 b_step_mask = _mm512_cmple_epu32_mask(b_vec.zmm, a_last_broadcasted);
                a += 16 - __lzcnt16((uint16_t)a_step_mask);
                b += 16 - __lzcnt16((uint16_t)b_step_mask);
            }

            // Handle the tail:
            c += intersect_u32_normal(a, b, a_end - a, b_end - b, out + c);
            // result += c; // And merge it with the main body result
            return c;
        }

        long double avg_skip = 0;
        long double avg_query_total_size = 0;
        long double avg_query_real_size = 0;

        // uint64_t total_query_duration = 0;

        void remove_duplicate_terms(std::vector<uint32_t>& terms) {
            std::sort(terms.begin(), terms.end());
            terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        }

        void query_test_intersection(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_round = 0;
            uint64_t avg_time_per = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    std::vector<la_vector_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const la_vector_enumerator <K> &a, const la_vector_enumerator <K> &b) {return a.n < b.n;});

                    int branch_num = query.size() == 2 ? 2 : 3;

                    for (int i = 0; i < branch_num; i++) {
                        index_sequences[i].query_init(decode_type, "intersection");
                    }

                    uint32_t query_id_idx = 0;
                    uint32_t candidate_posting = 0;
                    uint32_t equal_result = 0;
                    uint32_t candidate_posting_tmp = 0;
                    int intersection_size = query.size() == 2 ? index_sequences[1].n + 1 : index_sequences[2].n + 1;

                    K *intersection_result_p1, *intersection_result_p2;
                    std::vector<K, HugePageAllocator<K>> intersection_result_1(intersection_size); // for simd
                    std::vector<K, HugePageAllocator<K>> intersection_result_2(intersection_size);
                    // std::vector<K> intersection_result_1(intersection_size); // for normal
                    // std::vector<K> intersection_result_2(intersection_size);
                    intersection_result_p1 = intersection_result_1.data();
                    intersection_result_p2 = intersection_result_2.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p2[i] = i + 1;
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p1[i] = i + 1;
                    }

                    avg_time_per = 0;

                    if (query.size() == 1) {
                        throw std::runtime_error("query is empty");
                    }
                    else if (query.size() == 2) {
                        index_sequences[0].decode_query(intersection_result_p1, decode_type);
                        index_sequences[1].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n,  index_sequences[1].n, intersection_result_p1);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                    } else {
                        // first intersection
                        index_sequences[0].decode_query(intersection_result_p1, decode_type);
                        index_sequences[1].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, index_sequences[0].n, index_sequences[1].n, intersection_result_p1);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        // second intersection
                        index_sequences[2].decode_query(intersection_result_p2, decode_type);
                        avg_time_per += index_sequences[2].total_duration;

                        start = std::chrono::high_resolution_clock::now();
                        equal_result = intersect_u32_simd(intersection_result_p1, intersection_result_p2, equal_result, index_sequences[2].n, intersection_result_p1);
                        end = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        // next_geq candidate
                        if (index_sequences.size() > branch_num) {
                            for (query_id_idx = 3; query_id_idx < index_sequences.size(); query_id_idx++) {
                                index_sequences[query_id_idx].next_geq_init();
                                avg_time_per += index_sequences[query_id_idx].total_duration;
                            }

                            uint32_t intersection_length = equal_result;
                            equal_result = 0;
                            uint32_t intersection_idx = 0;
                            bool found = true;

                            start = std::chrono::high_resolution_clock::now();
                            while (intersection_idx < intersection_length) {
                                candidate_posting = intersection_result_p1[intersection_idx];
                                for (query_id_idx = 3; query_id_idx < index_sequences.size(); query_id_idx++) {
                                    candidate_posting_tmp = index_sequences[query_id_idx].nextgeq(candidate_posting);
                                    if (candidate_posting_tmp != candidate_posting) {
                                        candidate_posting = candidate_posting_tmp;
                                        found = false;
                                        break;
                                    }
                                }
                                if (found) {
                                    intersection_result_p2[equal_result++] = candidate_posting;
                                    intersection_idx++;
                                } else {
                                    intersection_idx++;
                                    // intersection_idx = lower_bound(intersection_result_p1 + intersection_idx + 1, intersection_result_p1 + intersection_length, candidate_posting) - intersection_result_p1;
                                    while (intersection_idx < intersection_length && intersection_result_p1[intersection_idx] < candidate_posting)
                                        intersection_idx++;
                                    found = true;
                                }
                            }

                            end = std::chrono::high_resolution_clock::now();
                            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                            avg_time_per += duration.count();
                        }
                    }


                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    avg_time_round += avg_time_per;

                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }

            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0ULL);

            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }

        void query_test_union(const std::vector<std::vector<uint32_t>> &query_list, const std::string &decode_type) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 1;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    std::vector<la_vector_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const la_vector_enumerator <K> &a, const la_vector_enumerator <K> &b) {return a.n < b.n;});

                    uint64_t union_size = 0, union_list_size = 0;
                    for (int i = 0; i < query.size(); i++) {
                        index_sequences[i].query_init(decode_type, "union");
                        union_size += index_sequences[i].n;
                        union_list_size = index_sequences[i].n > union_list_size ? index_sequences[i].n : union_list_size;
                    }

                    uint64_t equal_result = 0;
                    K *union_result_p1, *union_result_p2, *union_result_p3;
                    std::vector<K, HugePageAllocator<K>>  union_result_1(union_list_size); // for simd
                    std::vector<K, HugePageAllocator<K>>  union_result_2(union_size);
                    std::vector<K, HugePageAllocator<K>>  union_result_3(union_size);
                    // std::vector<K> union_result_1(union_list_size); // for normal
                    // std::vector<K> union_result_2(union_size);
                    // std::vector<K> union_result_3(union_size);
                    union_result_p1 = union_result_1.data();
                    union_result_p2 = union_result_2.data();
                    union_result_p3 = union_result_3.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < union_size; i++)
                            union_result_p3[i] = 0;
                        for (auto i = 0;i < union_size; i++)
                            union_result_p2[i] = 0;
                        for (auto i = 0;i < union_list_size; i++)
                            union_result_p1[i] = 0;
                    }

                    avg_time_per = 0;

                    if (query.size() == 1) {
                        throw std::runtime_error("query is empty");
                    }
                    if (query.size() == 2) {
                        index_sequences[0].decode_query(union_result_p1, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        index_sequences[1].decode_query(union_result_p2, decode_type);
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();
                    } else {
                        index_sequences[0].decode_query(union_result_p1, decode_type);
                        avg_time_per += index_sequences[0].total_duration;
                        index_sequences[1].decode_query(union_result_p2, decode_type);
                        avg_time_per += index_sequences[1].total_duration;

                        auto start = std::chrono::high_resolution_clock::now();
                        equal_result = union_u32_simd(union_result_p1, union_result_p2, index_sequences[0].n, index_sequences[1].n, union_result_p3);
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                        avg_time_per += duration.count();

                        for (int i = 2; i < query.size(); i++) {
                            index_sequences[i].decode_query(union_result_p1, decode_type);
                            avg_time_per += index_sequences[i].total_duration;

                            start = std::chrono::high_resolution_clock::now();
                            equal_result = union_u32_simd(union_result_p3, union_result_p1, equal_result, index_sequences[i].n, union_result_p2);
                            end = std::chrono::high_resolution_clock::now();
                            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                            avg_time_per += duration.count();

                            K *union_result_tmp = union_result_p2;
                            union_result_p2 = union_result_p3;
                            union_result_p3 = union_result_tmp;
                            // end_union_result = std::set_union(start_union_result, start_union_result + equal_result, index_sequences[i].current_value_vector.data(), index_sequences[i].current_value_vector.data()+index_sequences[i].n, start_union_result);
                        }
                    }

                    if (repeat == 0) {
                        total += equal_result;

                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;
                    }

                    avg_time_round += avg_time_per;
                    if (decode_type == "simd") {
                        for (auto &enumerator : index_sequences) {
                            enumerator.free_memory();
                        }
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }
            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0LL);

            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }

        static void intersection_candidate_test(std::vector<la_vector_enumerator <K>> &index_sequences, K* intersection_result_p1, uint32_t &equal_result, uint32_t &query_id_idx, uint32_t &candidate_posting_tmp, const uint32_t m) {
            uint32_t candidate_posting = index_sequences[0].docid();
            K* intersection_start = intersection_result_p1;

            while (candidate_posting < INT_MAX) {
                for (; query_id_idx < m; ++query_id_idx) { // find the same docID for all words
                    candidate_posting_tmp = index_sequences[query_id_idx].nextgeq(candidate_posting); // for word i  get next docID greater or equal to candidate
                    if (candidate_posting_tmp != candidate_posting) { // if exist a docID different from candidate
                        candidate_posting = candidate_posting_tmp; // candidate is the new docID
                        query_id_idx= 0; // i restarts
                        break; // break the loop and read the new candidate
                    }
                }
                if (query_id_idx == m) { // if all words have the same docID
                    *intersection_result_p1++ = candidate_posting; // add to the intersection
                    candidate_posting = index_sequences[0].nextgeq(candidate_posting + 1); // get a new candidate
                    query_id_idx = 1; // i restarts
                }

            }

            equal_result = intersection_result_p1 - intersection_start;
        }

        void query_test_intersection_benchmark(const std::vector<std::vector<uint32_t>> &query_list) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 5;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                int32_t query_id = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    query_id++;
                    if (query.size() < 2)
                        continue;

                    std::vector<la_vector_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const la_vector_enumerator <K> &a, const la_vector_enumerator <K> &b) {return a.n < b.n;});

                    uint32_t query_id_idx = 0;
                    uint32_t equal_result = 0;
                    uint32_t candidate_posting_tmp = 0;
                    const uint32_t m = index_sequences.size();
                    int32_t intersection_size = query.size() == 2 ? index_sequences[1].n + 1 : index_sequences[2].n + 1;

                    K *intersection_result_p1;
                    // std::vector<K, HugePageAllocator<K>> intersection_result_1(intersection_size); // for simd
                    // std::vector<K, HugePageAllocator<K>> intersection_result_2(intersection_size);
                    std::vector<K> intersection_result_1(intersection_size); // for normal
                    // std::vector<K> intersection_result_2(intersection_size);
                    intersection_result_p1 = intersection_result_1.data();
                    // intersection_result_p2 = intersection_result_2.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < intersection_size; i++)
                            intersection_result_p1[i] = i + 1 + k;
                    }

                    avg_time_per = 0;

                    for (int32_t i = 0; i < index_sequences.size(); i++) {
                        index_sequences[i].next_geq_init();
                        avg_time_per += index_sequences[i].total_duration;
                    }

                    auto start = std::chrono::high_resolution_clock::now();

                    intersection_candidate_test(index_sequences, intersection_result_p1, equal_result, query_id_idx, candidate_posting_tmp, m);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    avg_time_per += duration.count();

                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;

                    }

                    avg_time_round += avg_time_per;

                    for (auto &enumerator : index_sequences) {
                        enumerator.free_memory("normal");
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }
            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0LL);
            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }

        static void union_candidate_test(std::vector<la_vector_enumerator <K>> &index_sequences, K* result_p1, uint32_t &result, const uint32_t m) {
            uint32_t cur_doc = std::min_element(index_sequences.begin(), index_sequences.end(), [](la_vector_enumerator <K> lhs, la_vector_enumerator <K> rhs) {return lhs.docid() < rhs.docid();})->docid();
            K* result_start = result_p1;

            while (cur_doc < INT_MAX) {
                *result_p1++ = cur_doc;
                uint32_t next_doc = INT_MAX;
                for (int32_t i = 0; i < m; i++) {
                    if (index_sequences[i].docid() == cur_doc) {
                        index_sequences[i].next();
                    }
                    if (index_sequences[i].docid() < next_doc) {
                        next_doc = index_sequences[i].docid();
                    }
                }
                cur_doc = next_doc;
            }

            result = result_p1 - result_start;
        }

        void query_test_union_benchmark(const std::vector<std::vector<uint32_t>> &query_list) {
            std::vector<uint64_t> querys_per_times;
            uint64_t avg_time_per = 0;
            uint64_t avg_time_round = 0;
            K repeat_num = 5;

            for (K repeat = 0; repeat < repeat_num + 1; repeat++) {
                uint64_t total = 0;
                int32_t query_id = 0;
                avg_time_round = 0;
                for (auto &query : query_list) {
                    query_id++;

                    if (query.size() < 2)
                        continue;

                    std::vector<la_vector_enumerator <K>> index_sequences = load_model(query);

                    std::sort(index_sequences.begin(), index_sequences.end(), [](const la_vector_enumerator <K> &a, const la_vector_enumerator <K> &b) {return a.n < b.n;});


                    uint32_t equal_result = 0;
                    const uint32_t m = index_sequences.size();

                    K *union_result_p1;
                    std::vector<K> union_result_1(universe_size); // for normal
                    union_result_p1 = union_result_1.data();

                    // warm up
                    for (auto k = 0; k < 5; k++) {
                        for (auto i = 0;i < universe_size; i++)
                            union_result_p1[i] = i + 1 + k;
                    }

                    avg_time_per = 0;

                    for (int32_t i = 0; i < index_sequences.size(); i++) {
                        index_sequences[i].next_init();
                        avg_time_per += index_sequences[i].total_duration;
                    }

                    auto start = std::chrono::high_resolution_clock::now();

                    union_candidate_test(index_sequences, union_result_p1, equal_result, m);

                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    avg_time_per += duration.count();

                    if (repeat == 0) {
                        total += equal_result;
                        long double skip_tmp = 0;
                        long double size_tmp = 0;
                        for (int i = 0; i < index_sequences.size(); i++) {
                            skip_tmp += index_sequences[i].total_skip;
                            size_tmp += index_sequences[i].n;
                        }
                        avg_skip += skip_tmp / size_tmp;
                        avg_query_total_size += size_tmp;
                        avg_query_real_size += size_tmp - skip_tmp;

                    }

                    avg_time_round += avg_time_per;

                    for (auto &enumerator : index_sequences) {
                        enumerator.free_memory("simd");
                    }
                }
                if (repeat == 0)
                    std::cerr << "Total size: " << total << std::endl;
                if (repeat > 0)
                    querys_per_times.push_back(avg_time_round);
            }
            std::sort(querys_per_times.begin(), querys_per_times.end());
            uint64_t avg_time = std::accumulate(querys_per_times.begin(), querys_per_times.end(), 0LL);
            std::cerr << "Average query time: " <<  static_cast<long double> (avg_time) / query_list.size() / repeat_num / 1000 << ", " << "Median query time: " << static_cast<long double> (querys_per_times[querys_per_times.size() > 1 ? ((querys_per_times.size() + 1) / 2) : 0]) / query_list.size() / 1000 << std::endl;
            std::cerr << "Average skip rate: " << avg_skip / query_list.size() << ", Average query total size: " << avg_query_total_size / query_list.size() << ", Average query real size: " << avg_query_real_size / query_list.size() << std::endl;
        }


    public:
        K universe_size;
        void test_query(const std::string input_filename, const std::string &decode_type, const std::string &query_filename, const std::string query_type, const std::string &dataset_filename) {
            mm::file_source<K> input(dataset_filename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            universe_size = data[1];
            input.close();

            std::vector<std::vector<uint32_t>> query_list = read_query(query_filename);
            input_basename = input_filename;
            if (query_type == "AND")
                query_test_intersection(query_list, decode_type);
            else if (query_type == "OR")
                query_test_union(query_list, decode_type);
            else if (query_type == "ANDB")
                query_test_intersection_benchmark(query_list);
            else if (query_type == "ORB")
                query_test_union_benchmark(query_list);
            else {
                std::cerr << "Error: query_type must be AND, OR, ANDB or ORB" << std::endl;
            }
        }
    };
}

