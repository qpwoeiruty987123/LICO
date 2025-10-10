#pragma once
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <limits>

static inline long double LICO_S_scale = 200.0L;
static size_t global_m = 10;       // target max blocks for greedy
static size_t page_size = 1024;    // keys per page

static inline long double C_lower[10] = {1.97e8, 3.69e8, 5.78e8, 8.13e8, 1.09e9, 1.39e9, 1.76e9, 2.14e9, 2.86e9};
static inline long double C_value[10] = {0.152, 0.098, 0.065, 0.044, 0.031, 0.022, 0.017, 0.012, 0.009};

#ifndef MIN_PAGES_PER_BLOCK
#define MIN_PAGES_PER_BLOCK 1
#endif

using gap_prefix_type = __uint128_t;

// Half-open index ranges everywhere: [start_idx, end_idx)
struct Page {
    size_t start_idx = 0; // inclusive key index
    size_t end_idx   = 0; // exclusive key index
};

struct Block {
    size_t start_idx = 0; // inclusive key index
    size_t end_idx   = 0; // exclusive key index
    size_t epsilon   = 1; // epsilon*
    long double cost = 0.0L; // bits at epsilon*
};

struct BlockStats {
    size_t n_keys = 0;  // t - s
    size_t n_gaps = 0;  // n_keys - 1
    long double mu  = 0.0L; // mean gap
    long double var = 0.0L; // unbiased gap variance
};

static inline long double safe_log2l(long double x) {
    if (!(x > 0.0L)) x = LDBL_MIN; // guard
    // log2l is preferred, but log2 works fine with long double on many libcs
    return std::log2(x);
}

class DataPartition {
public:
    using K = uint32_t;

    long double LICO_C_coef = 0.09867L;

    // Inputs / state
    size_t total_size = 0;  // number of keys
    std::vector<gap_prefix_type> gap_prefix_sum;     // size n, prefix on gaps
    std::vector<gap_prefix_type> gap_prefix_squares; // size n, prefix on gap^2

    // Outputs
    std::vector<Page> pages;   // paging over keys, half-open
    std::vector<Block> blocks; // result blocks in key indices [s,t)
    std::string partition_type = "optimal";
    long double optimal_cost = 0.0L;
    long double greedy_cost  = 0.0L;

    DataPartition() = default;

    explicit DataPartition(const std::vector<K>& data) {
        total_size = data.size();
        build_page(data);
        build_prefix(data);
        calculate_C();
    }

    ~DataPartition() {
        std::vector<Block>().swap(blocks);
    }

    void calculate_C(){
        BlockStats block_sta = stats(pages[0].start_idx,  pages.back().end_idx);
        LICO_C_coef = 0;
        long double n_sigma_square = block_sta.var * block_sta.n_keys;
        if (n_sigma_square <= C_lower[8])
            LICO_C_coef = 0.09867L;
        // for (int  i = 0; i < 9; i++){
        //     if (n_sigma_square <= C_lower[i]) {
        //         // LICO_C_coef = C_value[i];
        //         LICO_C_coef = 0.09867L;
        //         break;
        //     }
        // }
        // std::cerr << n_sigma_square << "\tc_value: " << LICO_C_coef << std::endl;
    }

    void memory_clean() {
        std::vector<gap_prefix_type>().swap(gap_prefix_sum);
        std::vector<gap_prefix_type>().swap(gap_prefix_squares);
        std::vector<Page>().swap(pages);
    }

    void gap_prefix_clean() {
        std::vector<gap_prefix_type>().swap(gap_prefix_sum);
        std::vector<gap_prefix_type>().swap(gap_prefix_squares);
    }

    BlockStats stats(size_t s, size_t t) const {
        BlockStats r{};
        if (t <= s) return r;

        const size_t n_keys = t - s;
        r.n_keys = n_keys;

        if (n_keys < 2) return r;

        // gaps are between consecutive keys: indices [s+1 .. t-1]
        const size_t g0 = s;     // prefix index to subtract
        const size_t g1 = t - 1; // last gap index to include in prefix

        const long double sum =
            (long double)(gap_prefix_sum[g1] - gap_prefix_sum[g0]);
        const long double sum_sq =
            (long double)(gap_prefix_squares[g1] - gap_prefix_squares[g0]);

        const size_t n_gaps = n_keys - 1;
        r.n_gaps = n_gaps;

        const long double mu = sum / (long double)n_gaps;
        long double var = 0.0L;
        if (n_gaps >= 2) {
            var = (sum_sq - (long double)n_gaps * mu * mu) / (long double)(n_gaps - 1);
            if (var < 0.0L) var = 0.0L; // numeric guard
        }

        r.mu = mu;
        r.var = var;
        return r;
    }

    // 获取 size_t 类型中最高位的位置（0-based），x 必须 > 0
    static inline int highest_bit_position(size_t x) {
        return (int)(sizeof(size_t) * CHAR_BIT - 1) - __builtin_clzl(x);
    }

    // 返回：将 x 的最高位以下的所有位都置为 1
    static inline size_t fill_all_bits_from_msb(size_t x) {
        if (x == 0) return 0;
        int k = highest_bit_position(x);
        return (1UL << (k + 1)) - 1;
    }

    size_t epsilon_star_from(size_t n_keys, long double sigma2) {
        if (n_keys < 2 || !(sigma2 > 0.0L)) return 1;
        const long double n = (long double)(n_keys - 1); // #gaps
        const long double factor =
            (2.0L * std::log(2.0L)) * LICO_C_coef * LICO_S_scale * sigma2 * (n / (n + 1.0L));
        if (!(factor > 0.0L)) return 1;

        long double eps = std::sqrt(factor);
        if (!std::isfinite(eps) || eps <= 1.0L)
            return 1;
        const long double eps_up = std::ceil(eps);
        if (eps_up > (long double)std::numeric_limits<size_t>::max())
            return std::numeric_limits<size_t>::max();
        size_t eps_size_t = (size_t) eps_up;

        // return eps_size_t;

        // if (eps_size_t > 1 && LICO_C_coef == 0.001)
            // return fill_all_bits_from_msb(eps_size_t * 2);
        if (eps_size_t > 1)
            return fill_all_bits_from_msb(eps_size_t);
        else
            return 1; // epsilon* must be at least 1
    }

    size_t epsilon_star_range(size_t s, size_t t) {
        const auto st = stats(s, t);
        return epsilon_star_from(st.n_keys, st.var);
    }

    // S(P|ε*) ≈ n_gaps * ( 1.957 + 0.5 * log2(C*S*sigma^2) )
    long double bits_cost_star_from(size_t n_keys, long double sigma2) {
        if (n_keys < 2 || !(sigma2 > 0.0L)) return 0.0L;
        const size_t n_gaps = n_keys - 1;
        long double term = 1.957L + 0.5L * safe_log2l(LICO_C_coef * LICO_S_scale * sigma2);
        if (term < 0.0L) term = 0.0L; // never negative bits
        return (long double)n_gaps * term;
    }

    long double bits_cost_range(size_t s, size_t t) {
        const auto st = stats(s, t);
        return bits_cost_star_from(st.n_keys, st.var);
    }

    // ---------- Optimal (DP) over pages ----------
    void optimal_partition() {
        partition_type = "optimal";
        optimal_cost = 0.0L;
        blocks.clear();

        const size_t n = pages.size();
        if (n == 0) return;

        // dp[j]: min bits for covering page indices [0, j)
        std::vector<long double> dp(n + 1, 1e300L);
        std::vector<int> prev(n + 1, -1);
        dp[0] = 0.0L;

        for (size_t j = MIN_PAGES_PER_BLOCK; j <= n; ++j) {
            long double best = 1e300L;
            int best_i = -1;
            const size_t i_max = j - MIN_PAGES_PER_BLOCK;
            for (size_t i = 0; i <= i_max; ++i) {
                if (dp[i] >= 1e300L) continue;

                const size_t s_key = pages[i].start_idx;
                const size_t t_key = pages[j - 1].end_idx; // exclusive end of page j-1
                const long double cb = bits_cost_range(s_key, t_key);

                const long double cand = dp[i] + cb;
                if (cand < best) { best = cand; best_i = (int)i; }
            }
            dp[j] = best;
            prev[j] = best_i;
        }

        // Reconstruct
        int j = (int)n;
        while (j > 0) {
            const int i = prev[j];
            if (i < 0) {
                // Should not happen if transitions are correct; fallback: one big block
                Block blk;
                blk.start_idx = pages[0].start_idx;
                blk.end_idx   = pages[j - 1].end_idx;
                blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
                blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
                optimal_cost += blk.cost;
                blocks.push_back(blk);
                break;
            }
            Block blk;
            blk.start_idx = pages[(size_t)i].start_idx;
            blk.end_idx   = pages[(size_t)j - 1].end_idx;
            blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
            blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
            optimal_cost += blk.cost;
            blocks.push_back(blk);
            j = i;
        }
        std::reverse(blocks.begin(), blocks.end());
    }

    // ---------- Greedy (top-down) over pages ----------
    // Splits the most expensive interval until we have up to global_m blocks
    void greedy_partition() {
        partition_type = "greedy";
        greedy_cost = 0.0L;
        blocks.clear();

        using Interval = std::tuple<long double, size_t, size_t>; // (cost, s_page, t_page) with [s,t)
        auto cmp = [](const Interval& a, const Interval& b){ return std::get<0>(a) < std::get<0>(b); }; // max-heap
        std::priority_queue<Interval, std::vector<Interval>, decltype(cmp)> Q(cmp);

        const size_t n = pages.size();
        if (n == 0) return;

        // initial whole range
        {
            const size_t s_key = pages[0].start_idx;
            const size_t t_key = pages[n - 1].end_idx; // exclusive
            const long double full = bits_cost_range(s_key, t_key);
            Q.emplace(full, 0, n);
        }

        // helper to finalize a page interval [s,t)
        auto finalize_block = [&](size_t s, size_t t) {
            Block blk;
            blk.start_idx = pages[s].start_idx;
            blk.end_idx   = pages[t - 1].end_idx; // exclusive
            blk.epsilon   = epsilon_star_range(blk.start_idx, blk.end_idx);
            blk.cost      = bits_cost_range(blk.start_idx, blk.end_idx);
            greedy_cost  += blk.cost;
            blocks.push_back(blk);
        };

        if (LICO_C_coef != 0) {
            while (Q.size() < global_m && !Q.empty()) {
                auto [cur_cost, s, t] = Q.top();
                Q.pop();

                // number of pages in interval
                const size_t pn = t - s;
                // To split into [s,k) and [k,t), both must have at least MIN_PAGES_PER_BLOCK pages
                const size_t k_begin = s + MIN_PAGES_PER_BLOCK;
                const size_t k_end   = t - MIN_PAGES_PER_BLOCK;

                if (pn < 2 * MIN_PAGES_PER_BLOCK || k_begin > k_end) {
                    finalize_block(s, t);
                    continue;
                }

                // Find best split k
                size_t best_k = k_begin;
                long double best_sum = 1e300L;
                for (size_t k = k_begin; k <= k_end; ++k) {
                    const long double left  = bits_cost_range(pages[s].start_idx,     pages[k - 1].end_idx);
                    const long double right = bits_cost_range(pages[k].start_idx,     pages[t - 1].end_idx);
                    const long double sum_lr = left + right;
                    if (sum_lr < best_sum) { best_sum = sum_lr; best_k = k; }
                }

                // Push children intervals back to the heap
                {
                    const long double left_cost =
                        bits_cost_range(pages[s].start_idx, pages[best_k - 1].end_idx);
                    Q.emplace(left_cost, s, best_k);
                }
                {
                    const long double right_cost =
                        bits_cost_range(pages[best_k].start_idx, pages[t - 1].end_idx);
                    Q.emplace(right_cost, best_k, t);
                }
            }
        } else {
            LICO_C_coef = 0.001L;
        }


        // Finalize remaining intervals
        while (!Q.empty()) {
            auto [_, s, t] = Q.top();
            Q.pop();
            finalize_block(s, t);
        }

        // Sort blocks by key start
        std::sort(blocks.begin(), blocks.end(),
                  [](const Block& a, const Block& b){ return a.start_idx < b.start_idx; });
    }

    // ---------- Reporting ----------
    void summarize() const {
        std::cerr << "Partition type: " << partition_type << "\n";
        std::cerr << "Total Blocks: " << blocks.size() << "\n";
        long double total_bits = 0.0L;
        for (const auto &blk : blocks) {
            total_bits += blk.cost;
            std::cerr << "Block: [" << blk.start_idx << ", " << blk.end_idx
                      << "), eps*: " << blk.epsilon
                      << ", bits: " << (double)blk.cost << "\n";
        }
        std::cerr << "Total Bits: " << (double)total_bits << "\n";
    }

    void calculate_max_min_page_gap_variance(std::ofstream &file_gap_sta){
        long double max_gap_var = 0.0L;
        long double min_gap_var = 1e300L;
        for (const auto &page: pages) {
            auto page_stats = stats(page.start_idx, page.end_idx);
            if (page_stats.n_gaps == 0) continue;
            if (page_stats.var > max_gap_var) max_gap_var = page_stats.var;
            if (page_stats.var < min_gap_var) min_gap_var = page_stats.var;
        }
        file_gap_sta << "Max Page Gap Variance:\t" << (double)max_gap_var
                     << "\tMin Page Gap Variance:\t" << (double)min_gap_var
                     << "\tPage Number:\t" << pages.size() << "\t";
    }

    void summarize_gap_variance_per_page(std::ofstream &file_gap_sta) {
        file_gap_sta << pages.size() << "\n";
        for (const auto &page: pages) {
            auto page_stats = stats(page.start_idx, page.end_idx);
            // if (page_stats.n_gaps == 0) continue;
            file_gap_sta << page_stats.var << "\n";
        }
    }

private:
    // Build prefix sums of gaps (strictly increasing keys required)
    void build_prefix(const std::vector<K>& data) {
        const size_t n = data.size();
        gap_prefix_sum.assign(n, 0);
        gap_prefix_squares.assign(n, 0);
        if (n == 0) return;

        for (size_t i = 1; i < n; ++i) {
            if (data[i] <= data[i - 1])
                throw std::invalid_argument("Data must be strictly increasing (positive gaps).");
            const uint64_t g = data[i] - data[i - 1];
            gap_prefix_sum[i]      = gap_prefix_sum[i - 1] + (gap_prefix_type)g;
            gap_prefix_squares[i]  = gap_prefix_squares[i - 1]
                                   + (gap_prefix_type)g * (gap_prefix_type)g;
        }
    }

    // Build half-open pages over key indices
    void build_page(const std::vector<K>& data) {
        const size_t n = data.size();
        pages.clear();
        size_t s = 0;
        while (s < n) {
            size_t e = s + page_size;
            if (e > n) e = n; // exclusive
            pages.push_back({s, e});
            s = e;
        }
    }
};
