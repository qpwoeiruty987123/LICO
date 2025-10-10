#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <lico_index.hpp>
#include <lico_index_enumerate.hpp>
#include <lico_partition.hpp>
#include <../external/mm_file/include/mm_file/mm_file.hpp>

namespace lico_sequence
{

    template <typename K> // K is uint32_t or uint64_t
    class lico_estimater {
    public:
        std::vector<std::vector<LICOVariant>> index_partition_sequences; // store the partitioned index sequences
        std::vector<LICOVariant> index_sequences;
        uint64_t data_size = 0;
        uint64_t data_unequal = 0;
        uint64_t optPFD_size = 0;

        void test_epsilon(std::string input_basename) {
            std::map<uint32_t, uint32_t> epsilon_stats;
            std::cerr << "Read File [Build]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            std::cerr << "Universe Size: " << data[1] << std::endl;

            // random select 3 lists
            uint32_t list_num = 0;
            for (size_t i = 2; i < input.size();) {
                uint64_t n = data[i];
                list_num++;
                i += n + 1;
            }
            std::srand(list_num);
            std::vector<uint32_t> random_index;
            for (K i = 0; i < 3; i++)
                random_index.push_back(rand() % list_num);

            uint32_t list_idx = 0;
            for (size_t i = 2; i < input.size();){
                uint64_t n = data[i];

                ++list_idx;
                if (std::find(random_index.begin(), random_index.end(), list_idx) != random_index.end()) {
                    std::cerr << "List " << list_idx << "\t" << n << "\t";
                    std::vector<K> sequence(data + i + 1, data + i + n + 1);
                    std::vector<K> gaps;
                    for (size_t idx = 1; idx < sequence.size(); ++idx) {
                        gaps.push_back(sequence[idx] - sequence[idx - 1]);
                    }

                    long double mean = std::accumulate(gaps.begin(), gaps.end(), 0.0) / gaps.size();

                    long double variance = 0.0;
                    for (const K& gap : gaps) {
                        double diff = gap - mean;
                        variance += diff * diff;
                    }
                    variance /= gaps.size();

                    std::cerr << variance << "\t";

                    long double epsilon_star = std::sqrt((2.0 * std::log(2.0) * 0.109967L * 136.0 * static_cast<long double> (n) * variance) / static_cast<long double>(n + 1));

                    std::cerr << epsilon_star << "\n";

                    std::vector<LICOVariant> build_index_sequences;
                    build_index_sequences.reserve(20);
                    build_index_sequences.emplace_back(lico::LICO<K, 3>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 7>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 15>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 31>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 63>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 127>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 255>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 511>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 1023>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 2047>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 4095>(sequence));

                    for (auto &each_index: build_index_sequences) {
                        std::visit([](auto&& index){
                            std::cerr << index.Epsilon_Data << "\t" << static_cast<long double> (index.ground_truth_build_size_in_bytes()) / index.n * 8.0 << "\n";
                        }, each_index);
                    }
                }
                i += n + 1;
            }
            input.close();

        }

        void test_epsilon_uniform() {
            struct GapSpec {
                uint32_t gap_min;
                uint32_t gap_max;
                uint64_t size;
                uint64_t seed;
                const char* desc;
            };

            const uint64_t N = 1'000'000ULL;
            std::array<GapSpec, 3> specs{{
                GapSpec{1, 4,    N, 1234567ULL,   "U[1, 4]"},
                GapSpec{5, 25,   N, 9876543ULL,   "U[5, 25]"},
                GapSpec{100, 300, N, 20241009ULL, "U[100, 300]"},
            }};

            uint32_t list_idx = 0;

            for (const auto& s : specs) {
                ++list_idx;

                std::mt19937_64 rng(s.seed);
                std::uniform_int_distribution<uint32_t> dist(s.gap_min, s.gap_max);

                std::vector<K> sequence;
                sequence.reserve(static_cast<size_t>(s.size));
                K current = static_cast<K>(0);
                sequence.push_back(current);
                for (uint64_t i = 1; i < s.size; ++i) {
                    current = static_cast<K>(current + static_cast<K>(dist(rng)));
                    sequence.push_back(current);
                }

                std::vector<K> gaps;
                gaps.reserve(sequence.size() - 1);
                for (size_t i = 1; i < sequence.size(); ++i) {
                    gaps.push_back(static_cast<K>(sequence[i] - sequence[i - 1]));
                }

                long double mean = std::accumulate(gaps.begin(), gaps.end(), 0.0L) /
                                   static_cast<long double>(gaps.size());

                long double variance = 0.0L;
                for (const K& g : gaps) {
                    long double diff = static_cast<long double>(g) - mean;
                    variance += diff * diff;
                }
                variance /= static_cast<long double>(gaps.size());

                const uint64_t n = static_cast<uint64_t>(sequence.size());
                long double epsilon_star = std::sqrt(
                    (2.0L * std::log(2.0L) * 0.109967L * 136.0L *
                     static_cast<long double>(n) * variance) /
                    static_cast<long double>(n + 1));

                std::cerr << "List " << list_idx << " (" << s.desc << ")"
                          << "\t" << n
                          << "\t" << variance
                          << "\t" << epsilon_star
                          << "\n";


                std::vector<LICOVariant> build_index_sequences;
                build_index_sequences.reserve(11);
                build_index_sequences.emplace_back(lico::LICO<K, 3>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 7>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 15>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 31>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 63>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 127>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 255>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 511>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 1023>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 2047>(sequence));
                build_index_sequences.emplace_back(lico::LICO<K, 4095>(sequence));

                for (auto& each_index : build_index_sequences) {
                    std::visit([](auto& index) {
                        long double bits_per_element =
                            static_cast<long double>(index.ground_truth_build_size_in_bytes()) /
                            static_cast<long double>(index.n) * 8.0L;
                        std::cerr << index.Epsilon_Data << "\t" << bits_per_element << "\n";
                    }, each_index);
                }
            }
        }

    };
}
