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

    template <typename K, size_t epsilon = 64> // K is uint32_t or uint64_t
    class lico_builder {
    public:
        std::vector<std::vector<LICOVariant>> index_partition_sequences; // store the partitioned index sequences
        std::vector<LICOVariant> index_sequences;
        uint64_t data_size = 0;
        uint64_t data_unequal = 0;
        uint64_t optPFD_size = 0;

        void build_model(std::string input_basename, const uint32_t lico_m) {
            std::map<uint32_t, uint32_t> epsilon_stats;
            std::cerr << "\nEpsilon: " << epsilon << std::endl;
            std::cerr << "Read File [Build]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            assert(data[0] == 1);
            std::cerr << "Universe Size: " << data[1] << std::endl;
            long double total_cost = 0.0L;
            int total_unparted_blocks = 0;

            auto start = std::chrono::high_resolution_clock::now();

            int list_idx = 0;
            for (size_t i = 2; i < input.size();){
                uint64_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);

                if (sequence.back() == data[1])
                    throw std::runtime_error("Max is universe");

                ++list_idx;
                if (list_idx % 100 == 0)
                    std::cerr << list_idx << " ";

                if (epsilon == 0) { // divide data and build each partitioned list
                    lico_partition* data_partition = new lico_partition(sequence, lico_m);
                    data_partition -> greedy_partition();
                    total_cost += data_partition -> greedy_cost;


                    // data_partition -> optimal_partition();
                    // total_cost += data_partition -> optimal_cost;

                    // data_partition -> summarize();

                    data_partition -> memory_clean();

                    if (data_partition -> blocks.size() == 1)
                        total_unparted_blocks++;

                    std::vector<LICOVariant> build_index_sequences;
                    for (const auto & block: data_partition -> blocks) {
                        size_t Epsilon_Data = block.epsilon;
                        LICOVariant variant_index;

                        if (Epsilon_Data < 3) Epsilon_Data = 3;
                        if (Epsilon_Data > 1048575) Epsilon_Data = 1048575;
                        switch (Epsilon_Data) {
                            case 1: variant_index = lico::LICO<K, 1>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 3: variant_index = lico::LICO<K, 3>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 7: variant_index = lico::LICO<K, 7>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 15: variant_index = lico::LICO<K, 15>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 31: variant_index = lico::LICO<K, 31>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 63: variant_index = lico::LICO<K, 63>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 127: variant_index = lico::LICO<K, 127>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 255: variant_index = lico::LICO<K, 255>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 511: variant_index = lico::LICO<K, 511>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 1023: variant_index = lico::LICO<K, 1023>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 2047: variant_index = lico::LICO<K, 2047>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 4095: variant_index = lico::LICO<K, 4095>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 8191: variant_index = lico::LICO<K, 8191>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 16383: variant_index = lico::LICO<K, 16383>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 32767: variant_index = lico::LICO<K, 32767>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 65535: variant_index = lico::LICO<K, 65535>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 131071: variant_index = lico::LICO<K, 131071>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 262143: variant_index = lico::LICO<K, 262143>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 524287: variant_index = lico::LICO<K, 524287>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            case 1048575: variant_index = lico::LICO<K, 1048575>(sequence.begin() + block.start_idx, sequence.begin() + block.end_idx); break;
                            default: throw std::runtime_error("Unsupported epsilon Value: " + std::to_string(Epsilon_Data)); break;
                        }
                        epsilon_stats[Epsilon_Data]++;
                        build_index_sequences.emplace_back(variant_index);
                    }
                    index_partition_sequences.emplace_back(build_index_sequences);

                }
                else if (epsilon == 1) {                     // adapt each list build
                    std::vector<LICOVariant> build_index_sequences;
                    build_index_sequences.reserve(20);
                    build_index_sequences.emplace_back(lico::LICO<K, 1>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 3>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 7>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 15>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 31>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 63>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 127>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 255>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 511>(sequence));
                    build_index_sequences.emplace_back(lico::LICO<K, 1023>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 2047>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 4095>(sequence));

                    // build_index_sequences.emplace_back(lico::LICO<K, 8191>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 16383>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 32767>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 65535>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 131071>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 262143>(sequence));
                    // build_index_sequences.emplace_back(lico::LICO<K, 524287>(sequence));
                    uint64_t min_size = INT64_MAX;
                    uint32_t min_index = 0;
                    uint32_t min_epsilon = 0;
                    uint32_t index_count = 0;
                    for (auto &variant_index : build_index_sequences) {
                        std::visit([&min_size, &min_index, &index_count, &min_epsilon](auto &index) {
                            if (index.ground_truth_build_size_in_bytes() < min_size) {
                                min_size = index.ground_truth_build_size_in_bytes();
                                min_index = index_count;
                                min_epsilon = index.Epsilon_Data;
                            }
                        }, variant_index);
                        index_count++;
                    }
                    index_sequences.push_back(build_index_sequences[min_index]);
                    epsilon_stats[min_epsilon]++;
                    build_index_sequences.clear();
                }
                else {
                    index_sequences.push_back(lico::LICO<K, epsilon>(sequence));
                }
                data_size += sequence.size();
                i += n + 1;
            }
            input.close();

            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed= std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cerr << "\nBuild Time[ " << epsilon << " ]: " << elapsed.count() << " ms" << "\tTotal Cost: " << total_cost << "\tTotal Unparted Blocks: " << total_unparted_blocks<< std::endl;

            if (epsilon == 0) {
                int total_partitions = 0;
                for (auto it = epsilon_stats.begin(); it != epsilon_stats.end(); it++) {
                    total_partitions += it->second;
                    std::cerr << "Epsilon:\t" << it -> first << "\tCount:\t" << it -> second << std::endl;
                }
                std::cerr << "Total Partitions:\t" << total_partitions << " Epsilons Bytes: " << (total_partitions) << std::endl;
            } else if (epsilon == 1) {
                for (auto it = epsilon_stats.begin(); it != epsilon_stats.end(); it++) {
                    std::cerr << "Epsilon:\t" << it -> first  << "\tCount:\t" << it -> second << std::endl;
                }
            }
        }

        void statistic_index(std::string output_basename="", bool is_la_vector=false) {
            double segments_count = 0;
            uint64_t segments_count_real = 0;
            uint64_t segments_size = 0;
            uint64_t corrections_size = 0;
            uint64_t signs_size = 0;
            uint64_t errorpoint_size = 0;
            uint64_t total_list = 0;
            uint64_t calculate_by_total_size = 0;
            uint64_t slope_significand_max = 0;
            uint32_t slope_exponent_max = 0;
            double avg_covered = 0;

            std::ofstream file(output_basename + ".statistic_log.txt");
            std::ofstream file_r(output_basename + ".statistic_log_residual.txt");
            file << "The Epsilon with its segment counts under the list size of " << index_sequences.size() << std::endl;

            if (epsilon == 0) {
                K list_count = 0;
                total_list = index_partition_sequences.size();
                for (auto & index_partition : index_partition_sequences) {

                    uint64_t list_segments_size = 0, list_corrections_size = 0, list_signs_size = 0, list_size = 0;
                    for (auto& variant_index : index_partition) {
                        std::visit([&file_r, &list_count, &list_size, &avg_covered, &segments_count_real, &segments_size, &corrections_size, &signs_size, &errorpoint_size, &calculate_by_total_size, &slope_exponent_max, &slope_significand_max, &file, &list_segments_size, &list_corrections_size, &list_signs_size](auto &index) {
                            index.segment_init();
                            for (auto covered : index.seg_covered)
                                avg_covered += covered;
                            // slope_exponent_max = std::max(slope_exponent_max, index.segment_slope_exponent_max());
                            // slope_significand_max = std::max(slope_significand_max, index.segment_slope_significand_max());
                            segments_count_real += index.segments_size;

                            list_size += index.n;
                            list_segments_size += index.segment_size_in_bytes();
                            list_corrections_size += index.corrections_size_in_bytes();
                            list_signs_size += index.signs_size_in_bytes();

                            segments_size += index.segment_size_in_bytes();
                            corrections_size += index.corrections_size_in_bytes();
                            signs_size += index.signs_size_in_bytes();
                            errorpoint_size += index.errorPointCount;
                            file << "Epsilon:\t" << index.Epsilon_Data << "\tSegments Count:\t" << index.segments_size << "\tTotal bits per int:\t" << (segments_size + corrections_size + signs_size) * 8 / index.n << std::endl;
                            // index.normal_clean();
                            // index.normal_init();
                            // if (list_count == 12){
                                // index.report_residual_list(file_r);
                            // }
                        }, variant_index);
                    }
                    list_count++;
                    long double total_list_size = list_segments_size + list_corrections_size + list_signs_size;
                    file << "List:\t" << list_count << "\tPartition Number:\t" << index_partition.size() << "\tSegment Size:\t" << list_segments_size << "\tCorrection Size:\t" << list_corrections_size << "\tSign Size:\t" << list_signs_size << "\tList Length:\t" << list_size << "\tTotal bits per int:\t" << total_list_size / list_size * 8.0 << std::endl;
                }
                file_r.close();
            } else {
                total_list = index_sequences.size();
                K list_count = 0;
                for (auto& variant_index : index_sequences) {
                    list_count++;
                    uint64_t list_segments_size = 0, list_corrections_size = 0, list_signs_size = 0, list_size = 0;
                    std::visit([&output_basename, &is_la_vector, &list_count, &list_size, &avg_covered, &segments_count_real, &segments_size, &corrections_size, &signs_size, &errorpoint_size, &calculate_by_total_size, &slope_exponent_max, &slope_significand_max, &file, &list_segments_size, &list_corrections_size, &list_signs_size](auto &index) {
                    // index.normal_init();
                    index.segment_init();
                    for (auto &covered : index.seg_covered)
                        avg_covered += covered;
                    // slope_exponent_max = std::max(slope_exponent_max, index.segment_slope_exponent_max());
                    // slope_significand_max = std::max(slope_significand_max, index.segment_slope_significand_max());
                    segments_count_real += index.segments_size;
                    file << "Epsilon:\t" << index.Epsilon_Data << "\tSegments Count:\t" << index.segments_size;

                    list_size = index.n;
                    list_segments_size =  index.segment_size_in_bytes(is_la_vector);
                    list_corrections_size = index.corrections_size_in_bytes();
                    list_signs_size = index.signs_size_in_bytes();

                    segments_size += index.segment_size_in_bytes(is_la_vector);
                    corrections_size += index.corrections_size_in_bytes();
                    signs_size += index.signs_size_in_bytes();
                    errorpoint_size += index.errorPointCount;


                    if (list_count == 12 || list_count == 13 || list_count == 14 || list_count == 551 || list_count == 1061 || list_count == 2184 || list_count == 2270){
                        // std::cerr << list_count;
                        // index.normal_init();

                        // std::string file_name = output_basename + "_list[" + std::to_string(list_count) + "].bad_statistic_log_residual.txt";
                        // std::ofstream file_t(file_name);
                        // index.report_residual_list(file_t);
                        // file_t.close();
                    } else if (list_count == 648 || list_count == 964 || list_count == 489 ) {
                        // std::cerr << list_count;
                        // index.normal_init();

                        // std::string file_name = output_basename + "_list[" + std::to_string(list_count) + "].good_statistic_log_residual.txt";
                        // std::ofstream file_t(file_name);
                        // index.report_residual_list(file_t);
                        // file_t.close();
                    }

                }, variant_index);
                    long double total_list_size = list_segments_size + list_corrections_size + list_signs_size;
                    file << "List:\t" << list_count << "\tSegment Size:\t" << list_segments_size << "\tCorrection Size:\t" << list_corrections_size << "\tSign Size:\t" << list_signs_size << "\tList Length:\t" << list_size << "\tTotal bits per int:\t" << total_list_size / list_size * 8.0 << std::endl;
                }
            }
            segments_count = segments_count_real;
            double ratio = (segments_size + corrections_size + signs_size) / double(data_size) / double(sizeof(K));
            std::cerr << "Epsilon:\t" << epsilon << std::endl;
            std::cerr << "Integer Count:\t" << data_size << std::endl;
            std::cerr << "Segments count:\t" << segments_count_real << std::endl;
            std::cerr << "Error Point Count:\t" << errorpoint_size << std::endl;
            std::cerr << "Average Covered:\t" << avg_covered / segments_count << std::endl;
            std::cerr << "Average Segments per List:\t" << segments_count / total_list << std::endl;
            std::cerr << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
            std::cerr << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
            std::cerr << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
            std::cerr << "Compression Ratio:\t" << ratio << std::endl;
            long double total_size_in_bytes = segments_size + corrections_size + signs_size;
            long double total_size_in_gib = total_size_in_bytes / 1024.0 / 1024.0 / 1024.0;
            std::cerr << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int:" << std::endl;
            // std::cerr << "Slope Significand Max:\t" << slope_significand_max << " (" << BIT_WIDTH(slope_significand_max) << ")\tSlope Exponent Max:\t" << slope_exponent_max << " (" << BIT_WIDTH(slope_exponent_max) << ")" << std::endl;
            if (!output_basename.empty()){
                // std::ofstream file(output_basename + ".statistic_log.txt");
                file << "Epsilon:\t" << epsilon << std::endl;
                file << "Integer Count:\t" << data_size << std::endl;
                file << "Segments count:\t" << segments_count_real << std::endl;
                file << "Error Point Count:\t" << errorpoint_size << std::endl;
                // file << "Average Covered:\t" << avg_covered / segments_count << std::endl;
                file << "Average Length:\t" << segments_count / index_sequences.size() << std::endl;
                file << "Segment Size:\t" << segments_size << "\tbyte" << std::endl;
                file << "Corrections Size:\t" << corrections_size << "\tbyte" << std::endl;
                file << "Signs Size:\t" << signs_size << "\tbyte" << std::endl;
                file << "Compression Ratio:\t" << ratio << std::endl;
                file << "Total Size:\t" << total_size_in_bytes << "\tin bytes,\t" << total_size_in_gib << "\tin GiB,\t" << total_size_in_bytes / data_size << "\tbytes per int\t" << std::endl;
                // file << "Slope Significand Max:\t" << slope_significand_max << " (" << BIT_WIDTH(slope_significand_max) << ")\tSlope Exponent Max:\t" << slope_exponent_max << " (" << BIT_WIDTH(slope_exponent_max) << ")" << std::endl;
                std::cerr << "——Save statistic to: " << output_basename << ".statistic_log.txt——" << std::endl;
            }
        }

        void statistic_page_gap_variance_page_list(std::string input_basename, std::string output_basename) {
            std::ofstream file_gap_sta(output_basename + ".gap_variance_per_page.txt");
            std::cerr << "Read File [Page Gap]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();

            for (size_t i = 2; i < input.size();) {
                K n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);
                lico_partition* data_partition = new lico_partition(sequence);
                data_partition -> summarize_gap_variance_per_page(file_gap_sta);
                data_partition -> memory_clean();
                i += n + 1;
            }
            input.close();
        }

        void statistic_gap_list(std::string input_basename, std::string output_basename) {
            std::ofstream file_gap(output_basename + ".gaps.txt");
            std::ofstream file_gap_sta(output_basename + ".gaps_statistic.txt");
            std::cerr << "Read File [Gap]: " << input_basename << std::endl;
            mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::vector<double> gap_mean;
            std::vector<double> gap_variance;
            K idx = -1;
            bool save_gap = false;
            // random select 8 lists
            // std::vector<K> random_index;
            // for (K i = 0; i < 8; i++)
            //     random_index.push_back(rand() % input.size());
            for (size_t i = 2; i < input.size();) {
                K n = data[i];
                idx++;

                // if (std::find(random_index.begin(), random_index.end(), idx) != random_index.end())
                //     save_gap = true;
                // else
                //     save_gap = false;

                std::vector<K> sequence(data + i + 1, data + i + n + 1);

                if (save_gap)
                    file_gap << sequence.size() - 1 << std::endl;

                // statistic gap mean and variance
                K last = sequence[0];
                std::vector<int> gaps;
                uint64_t sum = 0;
                bool flag = true;
                for (int j = 1; j < sequence.size(); j++) {
                    // for (auto value : sequence) {
                    // if (flag) {
                    //     flag = false;
                    //     continue;
                    // }
                    K value = sequence[j];
                    int gap = value - last;

                    if (save_gap)
                        file_gap << gap << std::endl;

                    sum += gap;
                    gaps.push_back(gap);
                    last = value;
                }
                // std::cerr << n << " " <<  sequence.size() << " " << gaps.size() <<  std::endl;
                double mean = sum / double(gaps.size());
                // gap_mean.push_back(mean);
                double variance = 0;
                for (auto& gap : gaps) {
                    variance += (gap - mean) * (gap - mean);
                }
                variance /= gaps.size();
                // gap_variance.push_back(variance);

                lico_partition* data_partition = new lico_partition(sequence);
                data_partition -> calculate_max_min_page_gap_variance(file_gap_sta);
                data_partition -> memory_clean();
                file_gap_sta << "Key Min:\t" << sequence[0] << "\tKey Max:\t" << sequence.back() << "\t";
                file_gap_sta << "Gap Mean:\t" << mean << "\tGap Variance:\t" << variance << "\tCounts:\t" << n << std::endl;
                i += n + 1;
            }
            input.close();
            // double gap_mean_sum = 0;
            // double gap_variance_sum = 0;
            // for (auto& mean : gap_mean)
            // gap_mean_sum += mean;
            // for (auto& variance : gap_variance)
            // gap_variance_sum += variance;
            // file_gap_sta << "Total Mean:\t" << gap_mean_sum / gap_mean.size() << "\tTotal Variance:\t" << gap_variance_sum / gap_variance.size() << std::endl;
            // std::cerr << "Total Mean:\t" << gap_mean_sum / gap_mean.size() << "\tTotal Variance:\t" << gap_variance_sum / gap_variance.size() << std::endl;
            std::cerr << "Save Gaps List to: " << output_basename << ".gaps.txt" << std::endl;
            std::cerr << "Save Gaps Statistic to: " << output_basename << ".gaps_statistic.txt" << std::endl;
            file_gap.close();
            file_gap_sta.close();
        }

        void save_residual_and_residual_gaps_segments(std::string output_basename){
            std::ofstream file(output_basename + ".residual_and_residual_gaps_segments.txt");
            // random select 8 index
            std::vector<K> random_index;
            srand(42);
            for (K i = 0; i < 8; i++) {
                random_index.push_back(rand() % index_partition_sequences.size());
            }
            if (file.is_open()) {
                for (auto i : random_index) {
                    int block_idx = rand() % index_partition_sequences[i].size();
                    auto& variant_index = index_partition_sequences[i][block_idx];
                    std::visit([&i, &block_idx, &file](auto &index) {
                        index.normal_init();

                        int segment_idx = rand() % index.segments_size;
                        file << index.Epsilon_Data << "\t" << i << "\t" << block_idx << "\t" << segment_idx << "\n";

                        uint32_t first = 0, covered = index.seg_covered[segment_idx];
                        file << covered << "\n";
                        for (int j = 0; j < segment_idx; j++)
                            first += index.seg_covered[j];

                        int32_t last_residual = 0;
                        for (int j = first; j < first + covered; j++) {
                            int residual_gap = index.corrections_vector[j];
                            last_residual = residual_gap + last_residual;
                            uint32_t zigzag_residual_gap = zigzag(residual_gap);
                            file << last_residual << "\t" << residual_gap << "\t" << zigzag_residual_gap << "\n";
                        }

                    }, variant_index);
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void save_residual(const std::string output_basename) {
            std::ofstream file(output_basename + ".residual.txt");
            // random select 8 index
            std::vector<K> random_index;
            for (K i = 0; i < 8; i++) {
                random_index.push_back(rand() % index_sequences.size());
            }
            if (file.is_open()) {
                for (auto i : random_index) {
                    auto& variant_index = index_sequences[i];
                    std::visit([&file](auto &index) {
                        index.normal_init();
                        file << index.n << std::endl;
                        for (auto& correction : index.corrections_vector) {
                            file << correction << std::endl;
                        }
                    }, variant_index);
                }
                file.close();
            }
            else {
                std::cerr << "Couldn't open the log_file" << std::endl;
            }
        }

        void data_test(std::string input_basename) {
            K data_unequal_test = 0;
            if (epsilon == 0) {
                std::cerr << std::endl << "Data Test Epsilon: " << epsilon << std::endl;
                std::cerr << "Read File [Test]: " << input_basename << std::endl;
                mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
                K const* data = input.data();
                std::cerr << "Universe Size: " << data[1] << std::endl;
                K posi = 0;
                for (size_t i = 2; i < input.size();) {
                    K n = data[i];
                    std::vector<K> sequence(data + i + 1, data + i + n + 1);

                    lico_enumerator <K> enumerator_tmp = create_enumerator_from_indexes<K>(index_partition_sequences[posi++]);

                    std::vector<K> result_decode(enumerator_tmp.n);

                    // enumerator_tmp.normal_decode(result_decode.data());

                    enumerator_tmp.simd_init();
                    enumerator_tmp.simd_decode_512i(result_decode.data());

                    assert(result_decode.size() == sequence.size());

                    for (auto j = 0; j < result_decode.size(); j++) {
                        if (sequence[j] != result_decode[j]) {
                            data_unequal_test++;
                            // std::cerr << "Unequal Value: " << result_decode[j] << " " << sequence[j] << " " << posi << " " << j << std::endl;
                        }
                    }

                    i += n + 1;
                }

            }
            else {

                std::cerr << std::endl << "Data Test Epsilon: " << epsilon << std::endl;
                std::cerr << "Read File [Test]: " << input_basename << std::endl;
                mm::file_source<K> input(input_basename.c_str(), mm::advice::sequential);
                K const* data = input.data();
                std::cerr << "Universe Size: " << data[1] << std::endl;
                K posi = 0;
                for (size_t i = 2; i < input.size();) {
                    K n = data[i];
                    std::vector<K> sequence(data + i + 1, data + i + n + 1);
                    auto& variant_index = index_sequences[posi++];
                    std::visit([&sequence, &data_unequal_test](auto &index) {
                        index.segment_init();
                        std::vector<K> result1 = index.normal_decode();
                        assert(result1.size() == sequence.size());
                        for (auto j = 0; j < result1.size(); j++) {
                            if (sequence[j] != result1[j]) {
                                // std::cerr << sequence[j] << " " << result1[j] << " " << long(sequence[j]) - long(result1[j]) << std::endl;
                                data_unequal_test++;
                            }
                        }
                    }, variant_index);
                    i += n + 1;
                }
            }
            std::cerr << "Unequal postings: " << data_unequal_test << std::endl << std::endl;
        }

        void save_model(const std::string output_basename) {
            if (output_basename.empty() || output_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }
            std::cerr << "Save index to: " << output_basename << std::endl;

            // write_header
            std::ofstream out_header(output_basename + "idx.size", std::ios::binary);
            if (!out_header) {
                std::cerr << "Error: Cannot open idx.size for writing." << std::endl;
                return;
            }
            out_header.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));

            if (epsilon == 0) {
                K index_num = index_partition_sequences.size();
                out_header.write(reinterpret_cast<const char*>(&index_num), sizeof(index_num));
                // the size of each partition
                for (const auto& partition : index_partition_sequences) {
                    size_t partition_size = partition.size();
                    out_header.write(reinterpret_cast<const char*>(&partition_size), sizeof(size_t));
                }
                out_header.close();

                uint32_t index_count = 0;
                for (const auto& partition : index_partition_sequences) {
                    uint32_t partition_count = 0;
                    for (const auto& variant_index : partition) {
                        std::string filename = output_basename + std::to_string(index_count) + "_" + std::to_string(partition_count) + ".idx";
                        std::ofstream out(filename, std::ios::binary);
                        if (!out) {
                            std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
                            continue;
                        }
                        std::visit([&out](const auto& index) {
                            write_index_data(out, index);
                        }, variant_index);
                        out.close();
                        partition_count++;
                    }
                    index_count++;
                }
            } else {
                K index_num = index_sequences.size();
                out_header.write(reinterpret_cast<const char*>(&index_num), sizeof(index_num));
                out_header.close();

                K index_count = 0;
                for (const auto& variant_index : index_sequences) {
                    std::string filename = output_basename + std::to_string(index_count) + ".idx";
                    std::ofstream out(filename, std::ios::binary);
                    if (!out) {
                        std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
                        continue;
                    }
                    std::visit([&out](const auto& index) {
                        write_index_data(out, index);
                    }, variant_index);
                    out.close();
                    index_count++;
                }
            }
        }

        void load_model(const std::string input_basename) {
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

            in_header.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            if (epsilon == 0) {
                index_partition_sequences.clear();
                index_partition_sequences.reserve(index_num);

                for (K i = 0; i < index_num; ++i) {
                    size_t partition_size;
                    in_header.read(reinterpret_cast<char*>(&partition_size), sizeof(size_t));
                    std::vector<LICOVariant> partition;
                    partition.reserve(partition_size);

                    for (size_t j = 0; j < partition_size; ++j) {
                        std::string filename = input_basename + std::to_string(i) + "_" + std::to_string(j) + ".idx";
                        std::ifstream in(filename, std::ios::binary);
                        if (!in) {
                            std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                            continue;
                        }

                        uint32_t Epsilon_Data;
                        in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                        LICOVariant variant_index;
                        switch (Epsilon_Data) {
                            case 1: variant_index = lico::LICO<K, 1>(); break;
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
                            default:
                                std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl;
                                continue;
                        }

                        std::visit([&in, Epsilon_Data](auto& index) {
                            index.Epsilon_Data = Epsilon_Data;
                            read_index_data(in, index);
                        }, variant_index);

                        partition.emplace_back(std::move(variant_index));
                        in.close();
                    }
                    index_partition_sequences.emplace_back(std::move(partition));
                }
            } else {
                index_sequences.clear();
                index_sequences.reserve(index_num);

                for (K i = 0; i < index_num; ++i) {
                    std::string filename = input_basename + std::to_string(i) + ".idx";
                    std::ifstream in(filename, std::ios::binary);
                    if (!in) {
                        std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                        continue;
                    }

                    uint32_t Epsilon_Data;
                    in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));

                    LICOVariant variant_index;
                    switch (Epsilon_Data) {
                        case 1: variant_index = lico::LICO<K, 1>(); break;
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
                        default:
                            std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl;
                            continue;
                    }

                    std::visit([&in, Epsilon_Data](auto& index) {
                        index.Epsilon_Data = Epsilon_Data;
                        read_index_data(in, index);
                    }, variant_index);

                    index_sequences.push_back(std::move(variant_index));
                    in.close();
                }
            }
            in_header.close();
        }
    };
}
