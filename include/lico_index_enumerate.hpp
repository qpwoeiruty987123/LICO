#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <config.hpp>
#include <lico_index.hpp>
#include <tools.hpp>
#include <libdivide.h>

#if RESIDUAL_COMPRESS
#include <lico_fastpfor.hpp>
#endif

#ifndef HUGEPAGE
#define USE_HUGEPAGE 0
#else
#define USE_HUGEPAGE 1
#endif




namespace lico_sequence {
    using LICOIndex = lico::LICO<uint32_t>;

#if USE_HUGEPAGE
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
#endif

    template <typename K>
    class lico_enumerator{

        typedef int64_t Simd_Value;
        typedef uint64_t Slope_Value;
        typedef int32_t Correction_Value;
        typedef int32_t Segment_Value;

        public:

        struct segment{

            Segment_Value covered; // covered
            Segment_Value delta_y;
            Segment_Value delta_x;
            Segment_Value y_b; // 32 bits
            Segment_Value x_b;
            Segment_Value first;

            inline segment(Segment_Value first, Segment_Value delta_y, Segment_Value delta_x, Segment_Value y_b, Segment_Value x_b, Segment_Value covered) :
                    first(first), delta_y(delta_y), delta_x(delta_x), y_b(y_b), x_b(x_b), covered(covered) {}
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
        Segment_Value current_pos = 0;
        Segment_Value current_segment = 0;
        int32_t total_segment_size = 0;

#if USE_HUGEPAGE
        std::vector<K, HugePageAllocator<K>> current_value_vector;
#else
        std::vector<K> current_value_vector;
#endif


        void query_init(const std::string decode_type, const std::string query_type) {
            if (decode_type == "normal") {
                if (query_type == "intersection") {
                    current_pos = 0;
                    current_segment = 0;
                    total_skip = 0;
                } else if (query_type == "union") {
                    current_pos = 0;
                    current_segment = 0;
                    current_value = INT_MAX - 1;
                    total_skip = 0;
                }
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


        K* next_pointer;

        void next_init() {
            simd_init();
            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);

            current_pos = 0;
            current_value_vector.resize(n + 1);

            simd_decode_512i(current_value_vector.data());
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


        int64_t current_numerator;
        using branchfree_t = libdivide::branchfree_divider<int64_t>; // use libdivide to speed the division

        inline void build_segment_first_values() { // pre-decode the first value of each segment
            segment_first_values.clear();
            segment_first_values.resize(segments.size() + 1);

            auto start = std::chrono::high_resolution_clock::now();

            K* out = segment_first_values.data();
            size_t posi = 0;

            for (const auto& s : segments) {
                const int32_t j_x_diff = s.first - s.x_b;
                const int64_t delta_x_add = (s.delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                const int64_t numerator = static_cast<int64_t>(s.delta_y) * j_x_diff + delta_x_add;
        #if RESIDUAL_COMPRESS
                out[posi++] = static_cast<K>(numerator / s.delta_x + s.y_b + unzigzag_int32(corrections_vector[s.first]));
        #else
                out[posi++] = static_cast<K>(numerator / s.delta_x + s.y_b + corrections_vector[s.first]);
        #endif
            }
            out[posi] = INT_MAX; // sentinel

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
        }

        inline void next_geq_init() {
            residuals_decode();
            build_segment_first_values();

            total_segment_size = static_cast<int32_t>(segments.size());
            current_segment = 0;

            current_pos = 0;
            const segment* seg = &segments[0];
            current_correction = seg->y_b;

            const int32_t j_x_diff = seg->first - seg->x_b;
            const int64_t delta_x_add = (seg->delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
            current_numerator = static_cast<int64_t>(seg->delta_y) * j_x_diff + delta_x_add;

            current_value = segment_first_values[0];
            next_first_value = segment_first_values[1];
        }

        inline K nextgeq(K posting_value) {
            if (current_value >= posting_value) return current_value;

            const K* sfv = segment_first_values.data();

            // small sequential step
            if (posting_value >= next_first_value) {
                const int32_t seg_idx = current_segment;
                int32_t i = seg_idx + 1;
                const int32_t limit = std::min<int32_t>(seg_idx + 1 + 8, total_segment_size); // 8 times, if we still don't found the suitable segment, then use upper_bound

                while (i < limit && posting_value >= sfv[i]) {
                    ++i;
                }

                if (i < total_segment_size && posting_value < sfv[i]) {
                    current_segment = i - 1;
                } else {
                    const K* begin = sfv + i;
                    const K* end   = sfv + total_segment_size;
                    const K* it    = std::upper_bound(begin, end, posting_value);
                    current_segment = static_cast<int32_t>(it - sfv) - 1;

                    if (current_segment >= total_segment_size) {
                        current_pos = n;
                        current_value = INT_MAX;
                        return current_value;
                    }
                }

                const segment* seg = &segments[current_segment];
                current_pos = seg->first;
                current_correction = seg->y_b;
                const int32_t j_x_diff = seg->first - seg->x_b;
                const int64_t delta_x_add = (seg->delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                current_numerator = static_cast<int64_t>(seg->delta_y) * j_x_diff + delta_x_add;

                next_first_value = sfv[current_segment + 1];
            }

            // internal segment scan
            while (true) {
                const segment* seg = &segments[current_segment];
                const int32_t j_end   = seg->first + seg->covered;

                if (current_pos >= j_end) {
                    ++current_segment;
                    if (current_segment >= total_segment_size) {
                        current_pos = n;
                        // current_value = INT_MAX;
                        return INT_MAX;
                    }

                    const segment* nseg = &segments[current_segment];
                    current_pos = nseg->first;
                    current_correction = nseg->y_b;
                    const int32_t j_x_diff = nseg->first - nseg->x_b;
                    const int64_t delta_x_add = (nseg->delta_x >> 1) * ((j_x_diff > 0) - (j_x_diff < 0));
                    current_numerator = static_cast<int64_t>(nseg->delta_y) * j_x_diff + delta_x_add;
                    next_first_value = sfv[current_segment + 1];
                    continue;
                }

                const int32_t xb = seg->x_b;
                const int32_t xb_minus_1 = xb - 1;
                const int64_t dy = seg->delta_y;
                const int64_t dx_half64 = static_cast<int64_t>(seg->delta_x >> 1);
                const branchfree_t dx_div(seg->delta_x);

                const Correction_Value* corr_ptr = &corrections_vector[current_pos];
                const Correction_Value* corr_end = &corrections_vector[j_end];

                int32_t j = current_pos;

                {
                    const int32_t p1_end = std::min<int32_t>(j_end, xb_minus_1);
                    while (j < p1_end) {
        #if RESIDUAL_COMPRESS
                        current_correction += unzigzag_int32(*corr_ptr++);
        #else
                        current_correction += *corr_ptr++;
        #endif
                        const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                        current_numerator += dy;
                        if (v >= posting_value) {
                            current_value = v;
                            current_pos = ++j;
                            return current_value;
                        }
                        ++j;
                    }
                }

                // special process xb and xb - 1. Notably, here j is j + 1 for numerator, so we need check xb and xb - 1
                if ((j == xb_minus_1 || j == xb)) {
                    // if ((j == xb_minus_1 || j == xb) && j < j_end && corr_ptr < corr_end) {
                    for (int k = 0; k < 2  && (j == xb_minus_1 || j == xb); ++k) {
        #if RESIDUAL_COMPRESS
                        current_correction += unzigzag_int32(*corr_ptr++);
        #else
                        current_correction += *corr_ptr++;
        #endif
                        const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                        current_numerator += dy + dx_half64; // 两个特殊点都要额外加 dx_half
                        if (v >= posting_value) {
                            current_value = v;
                            current_pos = ++j;
                            return current_value;
                        }
                        ++j;
                    }
                }

                // last part
                while (j < j_end) {
        #if RESIDUAL_COMPRESS
                    current_correction += unzigzag_int32(*corr_ptr++);
        #else
                    current_correction += *corr_ptr++;
        #endif
                    const K v = static_cast<K>(current_numerator / dx_div + current_correction);
                    current_numerator += dy;
                    if (v >= posting_value) {
                        current_value = v;
                        current_pos = ++j;
                        return current_value;
                    }
                    ++j;
                }

                // go to next segment
                ++current_segment;
                if (current_segment >= total_segment_size) {
                    current_pos = n;
                    return INT_MAX;
                }

                const segment* nseg = &segments[current_segment];
                current_pos = nseg->first;
                current_correction = nseg->y_b;
                {
                    const int32_t j_x_diff2 = nseg->first - nseg->x_b;
                    const int64_t delta_x_add2 = (nseg->delta_x >> 1) * ((j_x_diff2 > 0) - (j_x_diff2 < 0));
                    current_numerator = static_cast<int64_t>(nseg->delta_y) * j_x_diff2 + delta_x_add2;
                }
                next_first_value = sfv[current_segment + 1];
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
            // uncompress residuals is not need to be decoded
        }


        // the scale decode algorithm
        void normal_decode(K* output) {
            K* __restrict out = output;

            Correction_Value* __restrict correct_pointer = corrections_vector.data();

            int32_t last_first = 0;

            const auto end_iter = segments.end();
            for (auto it = segments.begin(); it != end_iter; ++it) {
                const auto& seg = *it;

                const int32_t  covered = seg.covered;
                if (__builtin_expect(covered <= 0, 0)) { continue; }

                const int64_t  dy      = seg.delta_y;
                const int64_t  dx      = seg.delta_x;
                const int64_t  half    = dx >> 1;
                const int32_t  xb      = seg.x_b;
                const int32_t  yb      = seg.y_b;

                int64_t last_correction = static_cast<int64_t>(yb);
        #if RESIDUAL_COMPRESS
                last_correction += static_cast<int64_t>(unzigzag_int32(*correct_pointer++));
        #else
                last_correction += static_cast<int64_t>(*correct_pointer++);
        #endif

                int32_t j     = last_first;
                const int32_t j_end = j + covered;


                int64_t num = dy * static_cast<int64_t>(j - xb) + half * (static_cast<int64_t>(j > xb) - static_cast<int64_t>(j < xb));


                int64_t q = num / dx;
                int64_t r = num - q * dx;

                *out++ = static_cast<K>(q + last_correction);
                ++j;

                const int64_t step_q = dy / dx;
                const int64_t step_r = dy - step_q * dx;


                auto normalize_r = [&](int64_t& q_, int64_t& r_) {
                    if (__builtin_expect(r_ >= dx, 0)) {
                        r_ -= dx; ++q_;
                        if (r_ >= dx) { r_ -= dx; ++q_; }
                    } else if (__builtin_expect(r_ <= -dx, 0)) {
                        r_ += dx; --q_;
                        if (r_ <= -dx) { r_ += dx; --q_; }
                    }
                };


                auto sign_align = [&](int64_t& q_, int64_t& r_) {
                    const bool a_pos = (q_ > 0) || ((q_ == 0) && (r_ > 0));
                    const bool a_neg = (q_ < 0) || ((q_ == 0) && (r_ < 0));
                    if (__builtin_expect(a_pos && (r_ < 0), 0)) {
                        r_ += dx; --q_;
                    } else if (__builtin_expect(a_neg && (r_ > 0), 0)) {
                        r_ -= dx; ++q_;
                    }
                };

                // step without hafl
                auto step_nohalf = [&](){
        #if RESIDUAL_COMPRESS
                    last_correction += static_cast<int64_t>(unzigzag_int32(*correct_pointer++));
        #else
                    last_correction += static_cast<int64_t>(*correct_pointer++);
        #endif
                    q += step_q;
                    r += step_r;

                    normalize_r(q, r);
                    sign_align(q, r);

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                };

                // step with half
                auto step_half = [&](){
        #if RESIDUAL_COMPRESS
                    last_correction += static_cast<int64_t>(unzigzag_int32(*correct_pointer++));
        #else
                    last_correction += static_cast<int64_t>(*correct_pointer++);
        #endif
                    q += step_q;
                    r += step_r;
                    r += half;  // 只在这两次加 half

                    normalize_r(q, r);
                    sign_align(q, r);

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                };

                const int32_t e1 = (xb > j)     ? std::min(j_end, xb)     : j;      // [j, e1)
                const int32_t e2 = (xb + 1 > e1)? std::min(j_end, xb + 1) : e1;     // [e1+1, e2)

                // 1) [j, e1) no half
                while (j < e1) step_nohalf();

                // 2) j==xb, halft
                if (j < j_end && j == xb) step_half();

                // 3) (xb, e2) no half
                while (j < e2) step_nohalf();

                // 4) j==xb half
                if (j < j_end && j == xb + 1) step_half();

                // 5) last, no half
                while (j < j_end) step_nohalf();

                last_first = j_end;
            }
        }


        // used for SIMD
        constexpr static K align_val = 64; // for avx512
        std::vector<segment> segments_sort; // resorted segments
        std::vector<Segment_Value*> delta_y_simd;
        std::vector<Segment_Value*> delta_x_simd;
        std::vector<Segment_Value*> y_b_simd;
        std::vector<Segment_Value*> x_b_simd;
#if USE_HUGEPAGE
        std::vector<Correction_Value, HugePageAllocator<Correction_Value>> corrections_simd;
#else
        std::vector<Correction_Value> corrections_simd;
#endif
        std::vector<Correction_Value> corrections_tail;
        std::vector<Segment_Value*> first_simd;
        std::vector<Segment_Value*> covered_simd;
        std::vector<Segment_Value> cover_length;

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

#if USE_HUGEPAGE
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
#endif

        void free_memory(std::string decode_type = "simd") {
            if (decode_type == "simd") {
                for (auto i = 0; i < delta_y_simd.size(); i++) {
                    aligned_delete(delta_y_simd[i]);
                    aligned_delete(delta_x_simd[i]);
                    aligned_delete(y_b_simd[i]);
                    aligned_delete(x_b_simd[i]);
                    // aligned_delete(corrections_simd[i]);
                    // aligned_delete_huge(corrections_simd[i], cover_length[i]);
                    aligned_delete(first_simd[i]);
                    aligned_delete(covered_simd[i]);
                }
                std::vector<Segment_Value> ().swap(cover_length);


                std::vector<Correction_Value> ().swap(corrections_tail);
                std::vector<Segment_Value*> ().swap(delta_y_simd);
                std::vector<Segment_Value*> ().swap(delta_x_simd);
                std::vector<Segment_Value*> ().swap(y_b_simd);
                std::vector<Segment_Value*> ().swap(first_simd);
                std::vector<Segment_Value*> ().swap(covered_simd);
                std::vector<segment> ().swap(segments_sort);
#if USE_HUGEPAGE
                std::vector<K, HugePageAllocator<K>> ().swap(current_value_vector);
                std::vector<Correction_Value, HugePageAllocator<Correction_Value>> ().swap(corrections_simd);
#else
                std::vector<Correction_Value> ().swap(corrections_simd);
                std::vector<K> ().swap(current_value_vector);
#endif
            }

            std::vector<Correction_Value> ().swap(corrections_vector);
            std::vector<segment> ().swap(segments);
        }

        constexpr static  uint32_t simd_limit = 32;
        constexpr static double simd_group_limit = 0.0;

        void memory_layout() {
            std::vector<segment> tmp = segments;
            std::sort(tmp.begin(), tmp.end(), [](const segment& a, const segment& b) { return a.covered > b.covered;});

            std::vector<segment> simd;
            std::vector<segment> tail;
            simd.reserve(tmp.size());
            tail.reserve(tmp.size());

            const size_t n = tmp.size();
            size_t i = 0;

            while (i < n) {
                if (n - i < static_cast<size_t>(8)) {
                    tail.insert(tail.end(), tmp.begin() + static_cast<ptrdiff_t>(i), tmp.end());
                    break;
                }

                //  [i, i + 8)
                const size_t head = i;
                const size_t end  = i + static_cast<size_t>(8) - 1;

                const double head_cov = static_cast<double>(tmp[head].covered);
                const double tail_cov = static_cast<double>(tmp[end].covered);

                if (tail_cov >= head_cov * simd_group_limit && tail_cov >= simd_limit) {
                    simd.insert(simd.end(),
                                tmp.begin() + static_cast<ptrdiff_t>(i),
                                tmp.begin() + static_cast<ptrdiff_t>(i + 8));
                    i += static_cast<size_t>(8);
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


        // our Memory Re-Layout
        void simd_init() {
            idx = 0;
            memory_layout();

            Segment_Value min_cover = INT_MAX;
            Segment_Value max_min_covered = 0;

            uint64_t true_idx = 0;
            // for (auto it = segments_sort.begin(); it + 8 <= segments_sort.end(); it = it + 8) {
            assert(idx % 8 == 0);
            for (auto it = segments_sort.begin(); it +8 <= segments_sort.begin() + idx; it = it + 8) {
                alignas(align_val) Segment_Value *delta_y_simd_tmp = aligned_new<Segment_Value>(8);
                alignas(align_val) Segment_Value *delta_x_simd_tmp = aligned_new<Segment_Value>(8);
                alignas(align_val) Segment_Value *y_b_simd_tmp = aligned_new<Segment_Value>(8);
                alignas(align_val) Segment_Value *x_b_simd_tmp = aligned_new<Segment_Value>(8);
                alignas(align_val) Segment_Value *covered_tmp = aligned_new<Segment_Value>(8);
                alignas(align_val) Segment_Value *first_tmp = aligned_new<Segment_Value>(8);

                std::vector<segment> simd_segments(it, it + 8);
                true_idx += 8;

                int i = 0;
                for (auto its = simd_segments.begin(); its != simd_segments.end(); ++its, ++i) {
                    auto covered = its -> covered;
                    delta_y_simd_tmp[i] = static_cast<Segment_Value>(its -> delta_y);
                    delta_x_simd_tmp[i] = static_cast<Segment_Value>(its -> delta_x);
                    y_b_simd_tmp[i] = static_cast<Segment_Value>(its -> y_b);
                    x_b_simd_tmp[i] = static_cast<Segment_Value>(its -> x_b);
                    covered_tmp[i] = static_cast<Segment_Value> (covered);
                    first_tmp[i] = static_cast<Segment_Value> (its -> first);
                    min_cover = min_cover < covered ? min_cover : covered;
                }

                delta_y_simd.emplace_back(delta_y_simd_tmp);
                delta_x_simd.emplace_back(delta_x_simd_tmp);
                y_b_simd.emplace_back(y_b_simd_tmp);
                x_b_simd.emplace_back(x_b_simd_tmp);
                first_simd.emplace_back(first_tmp);
                covered_simd.emplace_back(covered_tmp);
                max_min_covered = min_cover;
                // max_min_covered = min_cover - min_cover % 2;
                cover_length.emplace_back(max_min_covered);
            }

            assert(idx == true_idx);

            residuals_decode();
            create_corrections();
            create_corrections_residual();
        }


        void create_corrections() {
            total_calculated = 0;
            for (int i = 0;i < cover_length.size(); i++) {
                total_calculated += cover_length[i] * 8;
            }
            corrections_simd.resize(total_calculated, 0);
            uint64_t corrections_pointer = 0;
            for (int i = 0;i < cover_length.size(); i++) {
                Segment_Value cover_length_tmp = cover_length[i];
                Segment_Value *first_tmp = first_simd[i];
                for (Segment_Value j = 0; j < cover_length_tmp; j++) {
                    for (Segment_Value k = 0; k < 8; k++) {
                        corrections_simd[corrections_pointer++] = corrections_vector[first_tmp[k] + j];
                    }
                }
            }
        }

        Segment_Value count_residual_size() {
            Segment_Value corrections_vector_residual_size = 0;
            for(int i = 0; i < cover_length.size(); i++) {
                Segment_Value *cover_tmp = covered_simd[i];
                Segment_Value cover_length_tmp = cover_length[i];
                for (Segment_Value k = 0; k < 8; k++) {
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
             corrections_tail.resize(count_residual_size());
             Segment_Value correct_pointers = 0;
             for(int i = 0; i < cover_length.size(); i++) {
                 Segment_Value *first_tmp = first_simd[i];
                 Segment_Value *cover_tmp = covered_simd[i];
                 Segment_Value cover_length_tmp = cover_length[i];
                 for (K k = 0; k < 8; k++) {
                     if (cover_length_tmp < cover_tmp[k]) {
                         for (Segment_Value pos = cover_length_tmp; pos < cover_tmp[k]; pos++){
                             corrections_tail[correct_pointers++] = corrections_vector[first_tmp[k] + pos];
                         }
                     }
                 }
             }

             for(auto it = segments_sort.begin() + idx; it < segments_sort.end(); it++) {
                 for (Segment_Value pos = it -> first; pos < (it -> first + it -> covered); pos++){
                     corrections_tail[correct_pointers++] = corrections_vector[pos];
                 }
             }
        }


        void simd_decode_512i(K* __restrict output) {
            Correction_Value* __restrict correct_pointer = corrections_tail.data();
            total_calculated = 0;
            total_calculated_add = 0;
            __m512i rerange_idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
            alignas(align_val) Correction_Value *vec_8x32int = aligned_new<Correction_Value>(8);
            alignas(align_val) Simd_Value *vec_8x64int = aligned_new<Simd_Value>(8);
            const Correction_Value *corrections_p = corrections_simd.data();
            const int cover_length_size = cover_length.size();

            auto start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < cover_length_size; i++) { // the slope is int64_t, the corrections and intercept are int32_t, we should align them
                Segment_Value  *first_tmp = first_simd[i];
                K *p0 = output + first_tmp[0];
                K *p1 = output + first_tmp[1];
                K *p2 = output + first_tmp[2];
                K *p3 = output + first_tmp[3];
                K *p4 = output + first_tmp[4];
                K *p5 = output + first_tmp[5];
                K *p6 = output + first_tmp[6];
                K *p7 = output + first_tmp[7];


                // #Block 1 predictions, we should calculate the first numerator
                // # 1 Load
                __m256i v_dx_32   = _mm256_load_epi32(delta_x_simd[i]);
                // cast to int64
                __m512i v_dy  = _mm512_cvtepi32_epi64(_mm256_load_epi32(delta_y_simd[i]));
                __m512d v_dx_pd  = _mm512_cvtepi32_pd(v_dx_32);
                __m512i v_dxd = _mm512_cvtepi32_epi64( _mm256_srli_epi32(v_dx_32, 1)); // delta_x / 2
                __m512i v_xb  = _mm512_cvtepi32_epi64(_mm256_load_epi32(x_b_simd[i]));
                __m512i v_j   = _mm512_cvtepi32_epi64(_mm256_load_epi32(first_simd[i]));

                //# 1 Calculate numerator
                // diff = j - x_b
                __m512i v_diff = _mm512_sub_epi64(v_j, v_xb);
                // (int64)delta_y * diff
                __m512i v_mul_dy = _mm512_mullo_epi64(v_dy, v_diff);

                // sgn(diff): 1 if diff>0, 0 if diff==0, -1 if diff<0
                __m512i v_zero_mask = _mm512_setzero_si512();
                __mmask8 m_pos = _mm512_cmp_epi64_mask(v_diff, v_zero_mask, _MM_CMPINT_GT); // diff > 0 mask
                __mmask8 m_neg = _mm512_cmp_epi64_mask(v_diff, v_zero_mask, _MM_CMPINT_LT); // diff < 0 mask
                __m512i v_sgn = _mm512_setzero_si512(); // dff = 0 mask
                v_sgn = _mm512_mask_set1_epi64(v_sgn, m_pos, 1); // set diff > 0
                v_sgn = _mm512_mask_set1_epi64(v_sgn, m_neg, -1); // set diff < 0, else = 0

                // delta_x_divide * sgn
                __m512i v_mul_dxd_sgn = _mm512_mullo_epi64(v_dxd, v_sgn);
                // the first numerator
                __m512i v_numerator = _mm512_add_epi64(v_mul_dy, v_mul_dxd_sgn);

                // #1 Calculate slope item
                __m512d v_numerator_pd = _mm512_cvtepi64_pd(v_numerator);
                v_numerator_pd = _mm512_div_pd(v_numerator_pd, v_dx_pd);
                __m256i v_int32 = _mm512_cvttpd_epi32(v_numerator_pd);


                // #1 residual block
                __m256i v_yb = _mm256_load_epi32(y_b_simd[i]); // y_b
                __m256i v_correction = _mm256_loadu_epi32(corrections_p); // residuals
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                __m256i mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals

                v_int32 = _mm256_add_epi32(v_int32, v_yb); // predictions + residuals + y_b

                // #1 store to lower 256 bits
                __m512i v_result = _mm512_castsi256_si512(v_int32); // lower 256 bits


                // #Block 2 predications
                // j = j +1
                __m512i v_one64 = _mm512_set1_epi64(1);
                v_j = _mm512_add_epi64(v_j, v_one64);

                // numerator + delta_y
                v_numerator = _mm512_add_epi64(v_numerator, v_dy);

                // (j == x_b or j == x_b+1)
                __mmask8 m_eq_xb = _mm512_cmpeq_epi64_mask(v_j, v_xb);
                __m512i v_xb_plus1 = _mm512_add_epi64(v_xb, v_one64);
                __mmask8 m_eq_xb_plus1 = _mm512_cmpeq_epi64_mask(v_j, v_xb_plus1);
                __mmask8 m_ind = m_eq_xb | m_eq_xb_plus1;

                // numerator + (delta_x_divide) * (j == x_b ? 1 : j == x_b + 1 ? 1 : 0)
                v_numerator= _mm512_mask_add_epi64(v_numerator, m_ind, v_numerator, v_dxd);

                // divid delta_x
                v_numerator_pd = _mm512_cvtepi64_pd(v_numerator);
                v_numerator_pd = _mm512_div_pd(v_numerator_pd, v_dx_pd);
                v_int32 = _mm512_cvttpd_epi32(v_numerator_pd);

                // #2 residual block
                v_correction = _mm256_loadu_epi32(corrections_p);
                corrections_p += 8;
#if RESIDUAL_COMPRESS
                mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                // predication + y_b + residuals
                v_int32 = _mm256_add_epi32(v_int32, v_yb);

                // store back
                v_result = _mm512_inserti64x4(v_result, v_int32, 1); // upper 256 bits
                v_result = _mm512_permutexvar_epi32(rerange_idx, v_result); // rearrange

                __m128i t0 = _mm512_extracti32x4_epi32(v_result, 0);
                __m128i t1 = _mm512_extracti32x4_epi32(v_result, 1);
                __m128i t2 = _mm512_extracti32x4_epi32(v_result, 2);
                __m128i t3 = _mm512_extracti32x4_epi32(v_result, 3);

                *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // b0, b1
                *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // c0, c1
                *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // d0, d1
                *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // e0, e1
                *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // f0, f1
                *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // g0, g1
                *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // h0, h1

                p0 = p0 + 2;
                p1 = p1 + 2;
                p2 = p2 + 2;
                p3 = p3 + 2;
                p4 = p4 + 2;
                p5 = p5 + 2;
                p6 = p6 + 2;
                p7 = p7 + 2;

                const Segment_Value cover_length_tmp = cover_length[i];
                total_calculated += cover_length_tmp; // test simd rate

                Segment_Value j_step = 2;
                const Segment_Value two_unit = cover_length_tmp - cover_length_tmp % 2;
                for (; j_step < two_unit; j_step += 2) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm512_add_epi64(v_j, v_one64);

                    // numerator + delta_y
                    v_numerator = _mm512_add_epi64(v_numerator, v_dy);

                    // (j == x_b or j == x_b+1)
                    m_eq_xb = _mm512_cmpeq_epi64_mask(v_j, v_xb);
                    m_eq_xb_plus1 = _mm512_cmpeq_epi64_mask(v_j, v_xb_plus1);
                    m_ind = m_eq_xb | m_eq_xb_plus1;

                    // numerator + (delta_x_divide) * (j == x_b ? 1 : j == x_b + 1 ? 1 : 0)
                    v_numerator= _mm512_mask_add_epi64(v_numerator, m_ind, v_numerator, v_dxd);

                    // divid delta_x
                    v_numerator_pd = _mm512_cvtepi64_pd(v_numerator);
                    v_numerator_pd = _mm512_div_pd(v_numerator_pd, v_dx_pd);
                    v_int32 = _mm512_cvttpd_epi32(v_numerator_pd);

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_int32, v_yb);

                    // tmp store to lower 256 bits
                    v_result = _mm512_castsi256_si512(v_int32);

                    // ——————————————————————————————————————————————
                    // #Decode Unit 2
                    v_j = _mm512_add_epi64(v_j, v_one64);

                    // numerator + delta_y
                    v_numerator = _mm512_add_epi64(v_numerator, v_dy);

                    // (j == x_b or j == x_b+1)
                    m_eq_xb = _mm512_cmpeq_epi64_mask(v_j, v_xb);
                    m_eq_xb_plus1 = _mm512_cmpeq_epi64_mask(v_j, v_xb_plus1);
                    m_ind = m_eq_xb | m_eq_xb_plus1;

                    // numerator + (delta_x_divide) * (j == x_b ? 1 : j == x_b + 1 ? 1 : 0)
                    v_numerator= _mm512_mask_add_epi64(v_numerator, m_ind, v_numerator, v_dxd);

                    // divid delta_x
                    v_numerator_pd = _mm512_cvtepi64_pd(v_numerator);
                    v_numerator_pd = _mm512_div_pd(v_numerator_pd, v_dx_pd);
                    v_int32 = _mm512_cvttpd_epi32(v_numerator_pd);

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_int32, v_yb);

                    // tmp store to upper 256bits
                    v_result = _mm512_inserti64x4(v_result, v_int32, 1); // upper 256 bits

                    // ——————————————————————————————————————————————
                    // #Store Unit 1
                    v_result = _mm512_permutexvar_epi32(rerange_idx, v_result); // rearrange

                    t0 = _mm512_extracti32x4_epi32(v_result, 0);
                    t1 = _mm512_extracti32x4_epi32(v_result, 1);
                    t2 = _mm512_extracti32x4_epi32(v_result, 2);
                    t3 = _mm512_extracti32x4_epi32(v_result, 3);


                    *((int64_t*)(p0)) = _mm_extract_epi64(t0, 0);  // a0, a1
                    *((int64_t*)(p1)) = _mm_extract_epi64(t0, 1);  // b0, b1
                    *((int64_t*)(p2)) = _mm_extract_epi64(t1, 0);  // c0, c1
                    *((int64_t*)(p3)) = _mm_extract_epi64(t1, 1);  // d0, d1
                    *((int64_t*)(p4)) = _mm_extract_epi64(t2, 0);  // e0, e1
                    *((int64_t*)(p5)) = _mm_extract_epi64(t2, 1);  // f0, f1
                    *((int64_t*)(p6)) = _mm_extract_epi64(t3, 0);  // g0, g1
                    *((int64_t*)(p7)) = _mm_extract_epi64(t3, 1);  // h0, h1

                    p0 = p0 + 2;
                    p1 = p1 + 2;
                    p2 = p2 + 2;
                    p3 = p3 + 2;
                    p4 = p4 + 2;
                    p5 = p5 + 2;
                    p6 = p6 + 2;
                    p7 = p7 + 2;

                }

                if (two_unit < cover_length_tmp) {
                    // #Decode Unit 1
                    // j = j +1
                    v_j = _mm512_add_epi64(v_j, v_one64);

                    // numerator + delta_y
                    v_numerator = _mm512_add_epi64(v_numerator, v_dy);

                    // (j == x_b or j == x_b+1)
                    m_eq_xb = _mm512_cmpeq_epi64_mask(v_j, v_xb);
                    m_eq_xb_plus1 = _mm512_cmpeq_epi64_mask(v_j, v_xb_plus1);
                    m_ind = m_eq_xb | m_eq_xb_plus1;

                    // numerator + (delta_x_divide) * (j == x_b ? 1 : j == x_b + 1 ? 1 : 0)
                    v_numerator= _mm512_mask_add_epi64(v_numerator, m_ind, v_numerator, v_dxd);

                    // divid delta_x
                    v_numerator_pd = _mm512_cvtepi64_pd(v_numerator);
                    v_numerator_pd = _mm512_div_pd(v_numerator_pd, v_dx_pd);
                    v_int32 = _mm512_cvttpd_epi32(v_numerator_pd);

                    // residual block
                    v_correction = _mm256_loadu_epi32(corrections_p);
                    corrections_p += 8;
#if RESIDUAL_COMPRESS
                    mask = _mm256_srai_epi32(_mm256_slli_epi32(v_correction, 31), 31); // unzigzag  mask = -(u & 1)
                    v_correction = _mm256_xor_si256(_mm256_srli_epi32(v_correction, 1), mask); // (u >> 1) ^ mask
#endif
                    v_yb = _mm256_add_epi32(v_yb, v_correction); // y_b + residuals, residuals is the gap

                    // predication + y_b + residuals
                    v_int32 = _mm256_add_epi32(v_int32, v_yb);

                    _mm256_store_epi32(vec_8x32int, v_int32);

                    *p0 = vec_8x32int[0];
                    *p1 = vec_8x32int[1];
                    *p2 = vec_8x32int[2];
                    *p3 = vec_8x32int[3];
                    *p4 = vec_8x32int[4];
                    *p5 = vec_8x32int[5];
                    *p6 = vec_8x32int[6];
                    *p7 = vec_8x32int[7];
                }

                // last residuals
                _mm256_store_epi32(vec_8x32int, v_yb);
                // last numerator
                _mm512_store_si512(vec_8x64int, v_numerator);

                const Segment_Value *cover_tmp = covered_simd[i];
                for (Segment_Value k = 0; k < 8; k++) {
                    const Segment_Value covered = cover_tmp[k];
                    if (cover_length_tmp < covered) {
                        K* output_tail = output + first_tmp[k] + cover_length_tmp;
                        int64_t numerator = vec_8x64int[k];
                        int32_t last_correction = vec_8x32int[k];
                        const int64_t delta_y = delta_y_simd[i][k];
                        const int32_t delta_x = delta_x_simd[i][k];
                        const int32_t delta_x_divide = delta_x >> 1; // delta_x / 2
                        const int32_t x_b = x_b_simd[i][k];

                        for (Segment_Value j = first_tmp[k] + cover_length_tmp; j < first_tmp[k] + covered; j++) {
#if RESIDUAL_COMPRESS
                            last_correction = last_correction + unzigzag_int32(*correct_pointer++);
#else
                            last_correction = last_correction + *correct_pointer++;
#endif
                            numerator += delta_y + (delta_x_divide) * ((j == x_b) | (j == x_b + 1));
                            *output_tail++ = numerator / delta_x + last_correction;
                        }
                    } else {
                        break;
                    }
                }
            }

            simd_tail_process(output, correct_pointer);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_duration += duration.count();
            total_calculated *= 8;
            aligned_delete(vec_8x32int);
            aligned_delete(vec_8x64int);
        }


        void simd_tail_process(K* __restrict output, Correction_Value* __restrict correct_pointer) {
            const auto end_iter = segments_sort.end();
            for (auto it = segments_sort.begin() + idx; it != end_iter; ++it) {
                const auto& seg = *it;

                const int32_t  covered = seg.covered;
                if (__builtin_expect(covered <= 0, 0)) { continue; }

                const int64_t  dy      = seg.delta_y;
                const int64_t  dx      = seg.delta_x;
                const int64_t  half    = dx >> 1;       // dx / 2
                const int32_t  xb      = seg.x_b;

                int32_t last_correction = (seg.y_b);
        #if RESIDUAL_COMPRESS
                last_correction += (unzigzag_int32(*correct_pointer++));
        #else
                last_correction += (*correct_pointer++);
        #endif

                int32_t j     = seg.first;
                K* __restrict out = output + j;
                const int32_t j_end = j + covered;

                // dy * (j - xb) + half * sign(j - xb)
                int64_t num = dy * (j - xb) + half * (static_cast<int64_t>(j > xb) - static_cast<int64_t>(j < xb));

                int64_t q = num / dx;
                int64_t r = num - q * dx;

                *out++ = static_cast<K>(q + last_correction);
                ++j;

                // num += dy
                const int64_t step_q = dy / dx;
                const int64_t step_r = dy - step_q * dx;

                auto normalize_r = [&](int64_t& q_, int64_t& r_) {
                    if (__builtin_expect(r_ >= dx, 0)) {
                        r_ -= dx; ++q_;
                        if (r_ >= dx) { r_ -= dx; ++q_; }
                    } else if (__builtin_expect(r_ <= -dx, 0)) {
                        r_ += dx; --q_;
                        if (r_ <= -dx) { r_ += dx; --q_; }
                    }
                };


                auto sign_align = [&](int64_t& q_, int64_t& r_) {
                    const bool a_pos = (q_ > 0) || ((q_ == 0) && (r_ > 0));
                    const bool a_neg = (q_ < 0) || ((q_ == 0) && (r_ < 0));
                    if (__builtin_expect(a_pos && (r_ < 0), 0)) {
                        r_ += dx; --q_;
                    } else if (__builtin_expect(a_neg && (r_ > 0), 0)) {
                        r_ -= dx; ++q_;
                    }
                };

                auto step_nohalf = [&](){
        #if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
        #else
                    last_correction += (*correct_pointer++);
        #endif
                    q += step_q;
                    r += step_r;

                    normalize_r(q, r);
                    sign_align(q, r);

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                };

                auto step_half = [&](){
        #if RESIDUAL_COMPRESS
                    last_correction += (unzigzag_int32(*correct_pointer++));
        #else
                    last_correction += (*correct_pointer++);
        #endif
                    q += step_q;
                    r += step_r;
                    r += half;

                    normalize_r(q, r);
                    sign_align(q, r);

                    *out++ = static_cast<K>(q + last_correction);
                    ++j;
                };

                const int32_t e1 = (xb > j) ? std::min(j_end, xb) : j;

                while (j < e1) step_nohalf();

                if (j < j_end && j == xb) step_half();

                if (j < j_end) step_half();

                while (j < j_end) step_nohalf();

            }
        }
    };


    template <typename K>
    static lico_enumerator<K> create_enumerator_from_indexes(std::vector<LICOIndex> &indexes) {
        lico_enumerator<K> enumerator;
        std::vector<int32_t> all_corrections_vector;
        std::vector<uint64_t> all_parted_size;
        std::vector<uint32_t> block_sizes_vector;
        std::vector<std::vector<uint32_t>> all_corrections_compress_fastpfor;

        uint64_t total_data_size = 0;
        int32_t last_first = 0;

        for (auto& idx : indexes) {
            idx.segment_init();
            total_data_size += idx.n;
            all_parted_size.push_back(idx.n);
            block_sizes_vector.push_back(idx.segments_size);

            for (int i = 0; i < idx.segments_size; ++i) {
                auto first = last_first;
                auto delta_y = idx.seg_delta_y[i];
                auto delta_x = idx.seg_delta_x[i];
                auto y_b = idx.seg_y_b[i];
                auto x_b = idx.seg_x_b[i];
                auto covered = idx.seg_covered[i];
                last_first += covered;
                enumerator.segments.emplace_back(first, delta_y, delta_x, y_b, x_b, covered);
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
    static lico_enumerator<K> create_enumerator_from_single_index(LICOIndex &idx) {
        lico_enumerator<K> enumerator;

        uint32_t last_first = 0;

        idx.segment_init();

        for (int i = 0; i < idx.segments_size; ++i) {
            auto first = last_first;
            auto delta_y = idx.seg_delta_y[i];
            auto delta_x = idx.seg_delta_x[i];
            auto y_b = idx.seg_y_b[i];
            auto x_b = idx.seg_x_b[i];
            auto covered = idx.seg_covered[i];
            last_first += covered;
            enumerator.segments.emplace_back(first, delta_y, delta_x, y_b, x_b, covered);
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

        return enumerator;
    }

}

