#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <config.hpp>
#include  <variant>
#include <lico_index.hpp>
#include <tools.hpp>

#if RESIDUAL_COMPRESS
#include <lico_fastpfor.hpp>
#endif


namespace lico_sequence {
    using LICOVariant = std::variant<
        lico::LICO<uint32_t, 0>,
        lico::LICO<uint32_t, 1>,
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

    template <typename T>
    class HugePageAllocator {
        constexpr static size_t PAGE_SIZE = 2 * 1024 * 1024; // 2MB
    public:
        using value_type = T;

        HugePageAllocator() = default;
        ~HugePageAllocator() = default;

        static T* allocate(std::size_t n) {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            // const size_t size = n * sizeof(T);
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            if (ptr == MAP_FAILED) {
                std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        static void deallocate(T* p, std::size_t n) noexcept {
            const size_t size = (n * sizeof(T) + PAGE_SIZE - 1) / PAGE_SIZE * PAGE_SIZE;
            munmap(p, size);
        }

        bool operator == (const HugePageAllocator&) const { return true; }
        bool operator != (const HugePageAllocator&) const { return false; }
    };

    template <typename K>
    class lico_enumerator{

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

        std::vector<std::vector<uint32_t>> corrections_compress_fastpfor;

        std::vector<uint32_t> block_sizes;

        void load_block_size(std::vector<uint32_t> all_block_sizes) {
            this -> block_sizes = std::move(all_block_sizes);
        }

        void load_residuals(uint64_t data_size, std::vector<Correction_Value> corrections_vector) {
            this -> n = data_size;
            // this -> corrections_vector = std::move(corrections_vector);
            this -> corrections_vector = std::move(corrections_vector);
        }

        void load_residuals_fastpfor(uint64_t data_size, std::vector<std::vector<uint32_t>> compress_fastpfor, std::vector<uint64_t> parted_sizes) {
            this -> n = data_size;
            this -> corrections_compress_fastpfor = std::move(compress_fastpfor);
            this -> parted_sizes = std::move(parted_sizes);
        }


        // used for Query Test
        K current_value = INT_MAX;
        K next_first_value = INT_MAX;
        Correction_Value current_correction = 0;
        Covered_Value current_pos = 0;
        Covered_Value current_segment = 0;
        uint32_t total_segment_size = 0;
        std::vector<K, HugePageAllocator<K>> current_value_vector;
        // std::vector<K> current_value_vector;

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
            } else if (decode_type == "simd") {
                simd_init();
                std::vector<Correction_Value> ().swap(corrections_vector);
                std::vector<segment> ().swap(segments);
                current_pos = 0;
                current_segment = 0;
                current_value = INT_MAX - 1;
                total_skip = 0;
            }
        }

        void decode_query(K* output, const std::string decode_type) {
            if (decode_type == "normal") {
                normal_decode(output);
            } else if (decode_type == "simd") {
                simd_decode_512i(output);
            }
        }

        long double total_skip = 0;

        void build_segment_first_values() {
            segment_first_values.clear();
            segment_first_values.resize(segments.size() + 1);
            K posi = 0;
            auto start = std::chrono::high_resolution_clock::now();
            for (auto &s : segments) {
#if RESIDUAL_COMPRESS
                segment_first_values[posi++] = s.intercept + unzigzag_int32(corrections_vector[s.first]);
#else
                segment_first_values[posi++] = s.intercept + corrections_vector[s.first];
#endif
            }
            segment_first_values[posi] = INT_MAX;
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        K* next_pointer;

        void next_init() {
            simd_init();
            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);

            current_pos = 0;
            current_value_vector.resize(n + 1);

            for (Covered_Value j = 0; j < n; j++) // warm up
                current_value_vector[j] = j + 1;

            simd_decode_512i(current_value_vector.data());
            current_value_vector[n] = INT_MAX;
            next_pointer = current_value_vector.data();
            current_value = *next_pointer;
        }

        K docid() {
            return current_value;
            // return *next_pointer;
        }

        void next() {
            // ++next_pointer;
            // if (current_pos < n - 1) {
            //     current_pos++;
            current_value = *++next_pointer;
            // }
            // else
                // current_value = INT_MAX;
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
                // next_first_value = segment_first_values[current_segment];
            }

            while (true){
                auto &seg = segments[current_segment];

                Covered_Value j_base = current_pos > seg.first ? (current_pos - seg.first) : 0;
                Slope_Value accum = (static_cast<Slope_Value>(j_base) * seg.slope_significand);
                const Correction_Value* corrections_pointer = corrections_vector.data() + j_base + seg.first;
                for (Covered_Value j = j_base; j < seg.covered; ++j) {
#if RESIDUAL_COMPRESS
                    current_correction += unzigzag_int32(*corrections_pointer++);
#else
                    current_correction += *corrections_pointer++;
#endif
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


        void residuals_decode() {
#if RESIDUAL_COMPRESS
            total_duration = 0;
            if (residual_compress_type == "fastpfor") {
                std::vector<uint32_t> uncompressed_output;
                uncompressed_output.resize(n);

                for (int32_t i = 0; i < 3; i++)
                    for (int32_t j = 0; j < n; j++)
                        uncompressed_output[j] = j + 1 + i; // just for warm up

                auto start = std::chrono::high_resolution_clock::now();

                decompress_residuals_fastpfor_parted(corrections_compress_fastpfor, uncompressed_output, parted_sizes);

                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                total_duration = duration.count();

                // actually, we can get int32 vector directly by the fastpfor decoding, beacuse 0 <= zigzag(residuals) <= 2*epsilon, and our epsilons are all smaller than (INT_MAX / 2), but dut to the FastPFor library only support ouput uint32_t type, so we cast them into int32_t
                corrections_vector.resize(n);
                for (int32_t i = 0; i < n; i++) {
                    corrections_vector[i] = static_cast<Correction_Value> (uncompressed_output[i]);
                }
            }
#endif
            // uncompress residuals is not need to decode
        }

        void normal_decode(K* output) {
            residuals_decode();
            Correction_Value* correction_pointer = corrections_vector.data();
            auto start = std::chrono::high_resolution_clock::now();
            const auto end_iter = segments.end();
            for (auto it = segments.begin(); it != end_iter; ++it) {
                const auto& seg = *it;
                Slope_Value significand =0;
                Correction_Value last_correction = seg.intercept;
                for (Covered_Value j = 0; j < seg.covered; ++j) {
#if RESIDUAL_COMPRESS
                    last_correction = last_correction + unzigzag_int32(*correction_pointer++);
#else
                    last_correction = last_correction + *correction_pointer++;
#endif
                    *output++ = static_cast<Correction_Value> (significand >> seg.slope_exponent) + last_correction;
                    significand += seg.slope_significand;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        // used for SIMD
        constexpr static K key_nums = 8;
        constexpr static K align_val = 64; // for avx512
        std::vector<segment> segments_sort; // resorted segments
        std::vector<Simd_Value*> slope_significand_simd;
        std::vector<Simd_Value*> slope_exponent_simd;
        std::vector<Intercept_Value*> intercept_simd;
        // std::vector<Correction_Value> corrections_simd;
        std::vector<Correction_Value, HugePageAllocator<Correction_Value>> corrections_simd;
        std::vector<Correction_Value> corrections_vector_residual;
        std::vector<Covered_Value*> first_simd;
        std::vector<Covered_Value*> covered_simd;
        std::vector<Covered_Value> cover_length;
        // std::vector<Covered_Value> tail_idx;
        // HashTable<Covered_Value, Covered_Value> decode_result_map;

        uint64_t total_calculated = 0;
        uint64_t total_calculated_add = 0;
        uint64_t conversion_time = 0;
        uint64_t total_duration = 0;

        uint64_t idx = 0;

        template <typename T>
        T* aligned_new(uint64_t num_elements) {
            void* ptr = std::aligned_alloc(align_val, num_elements * sizeof(T));
            if (!ptr) throw std::bad_alloc();
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete(T* ptr) {
            std::free(ptr);
        }

        constexpr static size_t HUGE_PAGE_SIZE = 2 * 1024 * 1024;

        template <typename T>
        T* aligned_new_huge(uint64_t num_elements) {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            void* ptr = mmap(nullptr, size,
                             PROT_WRITE | PROT_READ,
                             MAP_SHARED | MAP_ANONYMOUS | MAP_HUGETLB,
                             -1, 0);
            if (ptr == MAP_FAILED) {
                std::cerr << "Failed to allocate huge pages: " << std::strerror(errno) << " " << size << " " << num_elements << std::endl;
                throw std::bad_alloc();
            }
            return static_cast<T*>(ptr);
        }

        template <typename T>
        static void aligned_delete_huge(T* p, std::size_t num_elements) noexcept {
            const size_t size = (num_elements * sizeof(T) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE;
            munmap(p, size);
        }

        void free_memory(std::string decode_type = "simd") {
            if (decode_type == "simd") {
                for (auto i = 0; i < slope_significand_simd.size(); i++) {
                    aligned_delete(slope_significand_simd[i]);
                    aligned_delete(slope_exponent_simd[i]);
                    aligned_delete(intercept_simd[i]);
                    // aligned_delete(corrections_simd[i]);
                    // aligned_delete_huge(corrections_simd[i], cover_length[i]);
                    aligned_delete(first_simd[i]);
                    aligned_delete(covered_simd[i]);
                }
                std::vector<Covered_Value> ().swap(cover_length);
                // std::vector<Covered_Value> ().swap(tail_idx);
                std::vector<Correction_Value, HugePageAllocator<Correction_Value>> ().swap(corrections_simd);
                // std::vector<Correction_Value> ().swap(corrections_simd);
                std::vector<Correction_Value> ().swap(corrections_vector_residual);
                std::vector<Simd_Value*> ().swap(slope_significand_simd);
                std::vector<Simd_Value*> ().swap(slope_exponent_simd);
                std::vector<Intercept_Value*> ().swap(intercept_simd);
                std::vector<Covered_Value*> ().swap(first_simd);
                std::vector<Covered_Value*> ().swap(covered_simd);
                std::vector<segment> ().swap(segments_sort);
                std::vector<K, HugePageAllocator<K>> ().swap(current_value_vector);
                // aligned_delete(last_correction32);
                // vector<K> ().swap(current_value_vector);
            }

            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);
        }

        constexpr static  uint32_t simd_limit = 100;
        constexpr static double simd_group_limit = 0.9;

        void memory_layout() {
            // key_nums > 0
            assert(key_nums > 0);


            std::vector<segment> tmp = segments;
            std::sort(tmp.begin(), tmp.end(), [](const segment& a, const segment& b) { return a.covered > b.covered; });

            std::vector<segment> simd;
            std::vector<segment> tail;
            simd.reserve(tmp.size());
            tail.reserve(tmp.size());

            const size_t n = tmp.size();
            size_t i = 0;

            while (i < n) {
                if (n - i < static_cast<size_t>(key_nums)) {
                    tail.insert(tail.end(), tmp.begin() + static_cast<ptrdiff_t>(i), tmp.end());
                    break;
                }

                //  [i, i + key_nums)
                const size_t head = i;
                const size_t end  = i + static_cast<size_t>(key_nums) - 1;

                const double head_cov = static_cast<double>(tmp[head].covered);
                const double tail_cov = static_cast<double>(tmp[end].covered);

                if (tail_cov >= head_cov * simd_group_limit && tail_cov >= simd_limit) {
                    simd.insert(simd.end(),
                                tmp.begin() + static_cast<ptrdiff_t>(i),
                                tmp.begin() + static_cast<ptrdiff_t>(i + key_nums));
                    i += static_cast<size_t>(key_nums);
                } else {
                    tail.push_back(tmp[i]);
                    ++i;
                }
            }

            std::sort(tail.begin(), tail.end(), [](const segment& a, const segment& b) { return a.first < b.first; });

            segments_sort.clear();
            segments_sort.insert(segments_sort.end(), simd.begin(), simd.end());
            idx = simd.size();
            segments_sort.insert(segments_sort.end(), tail.begin(), tail.end());
        }


        // our SIMD
        void simd_init() {
            // segments_sort = segments;
            idx = 0;
            memory_layout();


            Covered_Value min_cover = INT_MAX;
            Covered_Value max_min_covered = 0;

            uint64_t true_idx = 0;
            // for (auto it = segments_sort.begin(); it + key_nums <= segments_sort.end(); it = it + key_nums) {
            assert(idx % key_nums == 0);
            for (auto it = segments_sort.begin(); it +key_nums <= segments_sort.begin() + idx; it = it + key_nums) {
                alignas(align_val) Simd_Value *slope_significand_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Simd_Value *slope_exponent_simd_tmp = aligned_new<Simd_Value>(key_nums);
                alignas(align_val) Intercept_Value *intercept_simd_tmp = aligned_new<Intercept_Value>(key_nums);
                alignas(align_val) Covered_Value *covered_tmp = aligned_new<Covered_Value>(key_nums);
                alignas(align_val) Covered_Value *first_tmp = aligned_new<Covered_Value>(key_nums);

                std::vector<segment> simd_segments(it, it + key_nums);
                // if (simd_segments.back().covered < simd_limit) {
                //     break;
                // }
                std::sort(simd_segments.begin(), simd_segments.end(), [](const segment &a, const segment &b) {return a.first < b.first;}); // part sorted
                true_idx += key_nums;

                int i = 0;
                for (auto its = simd_segments.begin(); its != simd_segments.end(); ++its, ++i) {
                    auto covered = its -> covered;
                    slope_significand_simd_tmp[i] = static_cast<Simd_Value>(its -> slope_significand);
                    slope_exponent_simd_tmp[i] = static_cast<Simd_Value>(its -> slope_exponent);
                    intercept_simd_tmp[i] = static_cast<Intercept_Value>(its -> intercept);
                    covered_tmp[i] = static_cast<Covered_Value> (covered);
                    first_tmp[i] = static_cast<Covered_Value> (its -> first);
                    min_cover = min_cover < covered ? min_cover : covered;
                }

                slope_significand_simd.emplace_back(slope_significand_simd_tmp);
                slope_exponent_simd.emplace_back(slope_exponent_simd_tmp);
                intercept_simd.emplace_back(intercept_simd_tmp);
                first_simd.emplace_back(first_tmp);
                covered_simd.emplace_back(covered_tmp);
                max_min_covered = min_cover - min_cover % 2;
                cover_length.emplace_back(max_min_covered);
            }

            // idx = true_idx;
            assert(idx == true_idx);
            // std::sort(segments_sort.begin() + idx, segments_sort.end(), [](const segment &a, const segment &b) {return a.first < b.first;});

            residuals_decode();
            create_corrections();
            create_corrections_residual();
        }


        void create_corrections() {
            total_calculated = 0;
            for (int i = 0;i < cover_length.size(); i++) {
                total_calculated += cover_length[i] * key_nums;
            }
            corrections_simd.resize(total_calculated, 0);
            uint64_t corrections_pointer = 0;
            for (int i = 0;i < cover_length.size(); i++) {
                Covered_Value cover_length_tmp = cover_length[i];
                Covered_Value *first_tmp = first_simd[i];
                for (Covered_Value j = 0; j < cover_length_tmp; j++) {
                    for (Covered_Value k = 0; k < key_nums; k++) {
                        corrections_simd[corrections_pointer++] = corrections_vector[first_tmp[k] + j];
                    }
                }
            }
        }

        Covered_Value count_residual_size() {
            Covered_Value corrections_vector_residual_size = 0;
            for(int i = 0; i < cover_length.size(); i++) {
                Covered_Value *cover_tmp = covered_simd[i];
                Covered_Value cover_length_tmp = cover_length[i];
                for (Covered_Value k = 0; k < key_nums; k++) {
                    auto covered = cover_tmp[k];
                    if (cover_length_tmp < covered)
                        corrections_vector_residual_size += covered - cover_length_tmp;
                }
            }
            for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++)
                corrections_vector_residual_size += it -> covered;
            return corrections_vector_residual_size;
        }

        void create_corrections_residual() {
             corrections_vector_residual.resize(count_residual_size());
             Covered_Value correct_pointers = 0;
             for(int i = 0; i < cover_length.size(); i++) {
                 Covered_Value *first_tmp = first_simd[i];
                 Covered_Value *cover_tmp = covered_simd[i];
                 Covered_Value cover_length_tmp = cover_length[i];
                 for (K k = 0; k < key_nums; k++) {
                     if (cover_length_tmp < cover_tmp[k]) {
                         for (Covered_Value pos = cover_length_tmp; pos < cover_tmp[k]; pos++){
                             corrections_vector_residual[correct_pointers++] = corrections_vector[first_tmp[k] + pos];
                         }
                     }
                 }
             }

             for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++) {
                 for (Covered_Value pos = it -> first; pos < (it -> first + it -> covered); pos++){
                     corrections_vector_residual[correct_pointers++] = corrections_vector[pos];
                 }
             }
        }

        

        void simd_decode_512i(K* output){
            Correction_Value* correct_pointer = corrections_vector_residual.data();
            total_calculated = 0;
            total_calculated_add = 0;
            __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            alignas(align_val) Correction_Value *last_correction32 = aligned_new<Correction_Value>(8);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = cover_length.size();

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them
                Covered_Value  *first_tmp = first_simd[i];
                K *p0 = output + first_tmp[0];
                K *p1 = output + first_tmp[1];
                K *p2 = output + first_tmp[2];
                K *p3 = output + first_tmp[3];
                K *p4 = output + first_tmp[4];
                K *p5 = output + first_tmp[5];
                K *p6 = output + first_tmp[6];
                K *p7 = output + first_tmp[7];

                // 0
                __m256i intercept_v = _mm256_load_epi32(intercept_simd[i]);

                __m256i corrections_v = _mm256_load_epi32(corrections_p);
#if RESIDUAL_COMPRESS
                __m256i mask = _mm256_srai_epi32(_mm256_slli_epi32(corrections_v, 31), 31); // unzigzag  mask = -(u & 1)
                corrections_v = _mm256_xor_si256(_mm256_srli_epi32(corrections_v, 1), mask); // (u >> 1) ^ mask
#endif
                corrections_p += 8;

                intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                __m512i result_v = _mm512_castsi256_si512(intercept_v); // lower 256 bits

                const Simd_Value *slope_significand_p = slope_significand_simd[i];
                const Simd_Value *slope_exponent_p = slope_exponent_simd[i];
                const __m512i slope_significand_v_tmp1 = _mm512_load_epi64(slope_significand_p);
                __m512i slope_significand_v1 = slope_significand_v_tmp1;
                __m512i slope_exponent_v1 = _mm512_load_epi64(slope_exponent_p);

                corrections_v = _mm256_load_epi32(corrections_p);
#if RESIDUAL_COMPRESS
                mask = _mm256_srai_epi32(_mm256_slli_epi32(corrections_v, 31), 31); // unzigzag  mask = -(u & 1)
                corrections_v = _mm256_xor_si256(_mm256_srli_epi32(corrections_v, 1), mask); // (u >> 1) ^ mask
#endif
                corrections_p += 8;

                intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                __m512i slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                __m256i int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                int32_v = _mm256_add_epi32(int32_v, intercept_v);
                result_v = _mm512_inserti64x4(result_v, int32_v, 1); // upper 256 bits
                result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                __m128i t0 = _mm512_extracti32x4_epi32(result_v, 0);
                __m128i t1 = _mm512_extracti32x4_epi32(result_v, 1);
                __m128i t2 = _mm512_extracti32x4_epi32(result_v, 2);
                __m128i t3 = _mm512_extracti32x4_epi32(result_v, 3);

                *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // a0, a1
                *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // a0, a1
                *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // a0, a1
                *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // a0, a1
                *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // a0, a1
                *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // a0, a1
                *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // a0, a1

                p0 = p0 + 2;
                p1 = p1 + 2;
                p2 = p2 + 2;
                p3 = p3 + 2;
                p4 = p4 + 2;
                p5 = p5 + 2;
                p6 = p6 + 2;
                p7 = p7 + 2;

                const Covered_Value cover_length_tmp = cover_length[i];
                total_calculated += cover_length_tmp;

                for (Covered_Value j = 2; j < cover_length_tmp; j += 2) {
                    corrections_v = _mm256_load_epi32(corrections_p);
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(corrections_v, 31), 31); // unzigzag  mask = -(u & 1)
                    corrections_v = _mm256_xor_si256(_mm256_srli_epi32(corrections_v, 1), mask); // (u >> 1) ^ mask
#endif
                    corrections_p += 8;

                    intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                    slope_significand_v1 = _mm512_add_epi64(slope_significand_v1, slope_significand_v_tmp1);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    int32_v = _mm256_add_epi32(int32_v, intercept_v);
                    result_v = _mm512_castsi256_si512(int32_v); // lower 256 bits

                    corrections_v = _mm256_load_epi32(corrections_p);
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(corrections_v, 31), 31); // unzigzag  mask = -(u & 1)
                    corrections_v = _mm256_xor_si256(_mm256_srli_epi32(corrections_v, 1), mask); // (u >> 1) ^ mask
#endif
                    corrections_p += 8;

                    intercept_v = _mm256_add_epi32(intercept_v, corrections_v);
                    slope_significand_v1 = _mm512_add_epi64(slope_significand_v1, slope_significand_v_tmp1);
                    slope_correct_v = _mm512_srlv_epi64(slope_significand_v1, slope_exponent_v1);
                    int32_v = _mm512_cvtepi64_epi32(slope_correct_v);
                    int32_v = _mm256_add_epi32(int32_v, intercept_v);
                    result_v = _mm512_inserti64x4(result_v, int32_v, 1); // upper 256 bits

                    result_v = _mm512_permutexvar_epi32(rerange_idx, result_v);

                    t0 = _mm512_extracti32x4_epi32(result_v, 0);
                    t1 = _mm512_extracti32x4_epi32(result_v, 1);
                    t2 = _mm512_extracti32x4_epi32(result_v, 2);
                    t3 = _mm512_extracti32x4_epi32(result_v, 3);

                    *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                    *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // a0, a1
                    *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // a0, a1
                    *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // a0, a1
                    *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // a0, a1
                    *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // a0, a1
                    *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // a0, a1
                    *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // a0, a1

                    p0 = p0 + 2;
                    p1 = p1 + 2;
                    p2 = p2 + 2;
                    p3 = p3 + 2;
                    p4 = p4 + 2;
                    p5 = p5 + 2;
                    p6 = p6 + 2;
                    p7 = p7 + 2;
                }

                _mm256_store_epi32(last_correction32, intercept_v);

                const Covered_Value *cover_tmp = covered_simd[i];
                for (uint8_t k = 0; k < 8; k++) {
                    const Covered_Value covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        const Slope_Value slope_significand_tmp = static_cast<Slope_Value> (slope_significand_p[k]);
                        Slope_Value slope_significand = slope_significand_tmp * cover_length_tmp;
                        const uint8_t slope_exponent = static_cast<uint8_t> (slope_exponent_p[k]);
                        Correction_Value last_correction = last_correction32[k];
                        K* p_tmp = output + first_tmp[k] + cover_length_tmp;
                        for (Covered_Value pos = cover_length_tmp; pos < covered; pos++) {
#if RESIDUAL_COMPRESS
                            last_correction = last_correction + unzigzag_int32(*correct_pointer++);
#else
                            last_correction = last_correction + *correct_pointer++;
#endif
                            *p_tmp++ = static_cast<Correction_Value> (slope_significand >> slope_exponent) + last_correction;
                            slope_significand += slope_significand_tmp;
                        }
                    }
                }
            }

            const auto end_iter = segments_sort.end();
            for (auto it = segments_sort.begin() + idx; it < end_iter; ++it) {
                const auto& seg = *it;
                K* dst = output + seg.first;
                Correction_Value last_correction = seg.intercept;
                Slope_Value slope_significand_tmp = 0;
                for (Covered_Value pos = 0; pos < seg.covered; ++pos) {
#if RESIDUAL_COMPRESS
                    last_correction = last_correction + unzigzag_int32(*correct_pointer++);
#else
                    last_correction = last_correction + *correct_pointer++;
#endif
                    *dst++ = static_cast<Correction_Value> (slope_significand_tmp >> seg.slope_exponent) + last_correction;
                    slope_significand_tmp += seg.slope_significand;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            total_calculated *= key_nums;
            aligned_delete(last_correction32);
        }
    
    };

    template <typename K>
    static lico_enumerator<K> create_enumerator_from_indexes(std::vector<LICOVariant> &indexes) {
        lico_enumerator<K> enumerator;
        std::vector<int32_t> all_corrections_vector;
        std::vector<uint64_t> all_parted_size;
        std::vector<uint32_t> block_sizes_vector;
        std::vector<std::vector<uint32_t>> all_corrections_compress_fastpfor;

        uint64_t total_data_size = 0;
        uint32_t last_first = 0;

        for (auto& index : indexes) {
            std::visit([&enumerator, &all_parted_size, &total_data_size, &last_first, &all_corrections_vector, &all_corrections_compress_fastpfor, &block_sizes_vector](auto &idx) {
                idx.segment_init();
                total_data_size += idx.n;
                all_parted_size.push_back(idx.n);
                block_sizes_vector.push_back(idx.segments_size);

                for (int i = 0; i < idx.segments_size; ++i) {
                    auto first = last_first;
                    auto intercept = idx.seg_intercept[i];
                    auto slope_exponent = idx.seg_slope_exponent[i];
                    auto slope_significand = idx.seg_slope_significand[i];
                    auto covered = idx.seg_covered[i];
                    last_first += covered;
                    enumerator.segments.emplace_back(first, intercept, slope_exponent, slope_significand, covered);
                }

#if RESIDUAL_COMPRESS
                if (residual_compress_type == "fastpfor") {
                    all_corrections_compress_fastpfor.push_back(idx.corrections_compress);
                } else {
                    throw std::runtime_error("residual_compress_type not recognised");
                }
#else
                    idx.residual_init();
                    all_corrections_vector.insert(all_corrections_vector.end(), idx.corrections_vector.begin(), idx.corrections_vector.end());
#endif
            }, index);
        }

#if RESIDUAL_COMPRESS
        if (residual_compress_type == "fastpfor") {
            enumerator.load_residuals_fastpfor(total_data_size, all_corrections_compress_fastpfor, all_parted_size);
        } else {
            throw std::runtime_error("residual_compress_type not recognised");
        }
#else
        enumerator.load_residuals(total_data_size, all_corrections_vector);
#endif

        enumerator.load_block_size(block_sizes_vector);

        return enumerator;
    }

    template <typename K>
    static lico_enumerator<K> create_enumerator_from_single_index(LICOVariant &index) {
        lico_enumerator<K> enumerator;

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

#if RESIDUAL_COMPRESS
            if (residual_compress_type == "fastpfor") {
                std::vector<uint64_t> parted_size{idx.n};
                std::vector<std::vector<uint32_t>> compress_fastpfor{idx.corrections_compress};
                enumerator.load_residuals_fastpfor(idx.n, compress_fastpfor, parted_size);
            } else {
                throw std::runtime_error("residual_compress_type not recognised");
            }
#else
                idx.residual_init();
                enumerator.load_residuals(idx.n, idx.corrections_vector);
#endif
            enumerator.load_block_size(std::vector<uint32_t>{idx.segments_size});
            idx.free_memory();
        }, index);

        return enumerator;
    }

}

