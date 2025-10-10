#pragma once

#include <iostream>
#include <cassert>
#include <vector>
#include <lico_index.hpp>
#include <lico_index_enumerate.hpp>
// #include <perf_event.hpp>

namespace lico_sequence {
    template <typename K, size_t epsilon = 64> // K is uint32_t or uint64_t
    class lico_decoder{

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

        void decode_test(lico_enumerator <K> &index, std::string &decode_type, bool warm_up = true) {
            total_list++;

            if (decode_type == "simd") {
                index.simd_init();

                // std::vector<K> result1(index.n);
                std::vector<K, HugePageAllocator<K>> result1(index.n);
                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result1[i] = 0;
                        }
                    }
                }
                // perf test
                // PerfEvent perf_event;
                // perf_event.startCounters();

                index.simd_decode_512i(result1.data());

                // perf_event.stopCounters();
                // perf_add(perf_event);

                // debug

                // std::vector<K> result2(index.n);
                // index.normal_decode(result2.data());
                //
                // for (int i = 0; i < result1.size(); i ++) {
                //     assert(result1[i] == result2[i]);
                // }
                // std::vector<K> ().swap(result2);


                // std::vector<K> ().swap(result1);
                std::vector<K, HugePageAllocator<K>> ().swap(result1);

                total_calculated += index.total_calculated;

                if (index.total_duration > max_decode_time)
                    max_decode_time = index.total_duration;
                if (index.total_duration < min_decode_time)
                    min_decode_time = index.total_duration;

                // std::cerr << " list id: " << total_list << " list size: " << index.n << " cost time: " << index.total_duration << " simd_rate: " << double(index.total_calculated) / index.n << " ns per int: " << double(index.total_duration) / index.n * 1000 << " segment number: " << index.segments.size() << std::endl;
            }
            else if (decode_type == "normal") {
                // std::vector<K> result2(index.n);
                std::vector<K, HugePageAllocator<K>> result2(index.n);
                if (warm_up) {
                    for (int warm_time = 0; warm_time < 5; warm_time++) {
                        for (int i = 0; i < index.n; i ++) {
                            result2[i] = 0;
                        }
                    }
                }
                // perf test
                // PerfEvent perf_event;
                // perf_event.startCounters();

                index.normal_decode(result2.data());

                // perf_event.stopCounters();
                // perf_add(perf_event);

                // std::vector<K> ().swap(result2);
                std::vector<K, HugePageAllocator<K>> ().swap(result2);

                if (index.total_duration > max_decode_time)
                    max_decode_time = index.total_duration;
                if (index.total_duration < min_decode_time)
                    min_decode_time = index.total_duration;
            }

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

            // std::cerr << "Performance: " << std::endl;
            // for (const auto &header : perf_header) {
            //     std::cerr << std::setw(15) << header;
            // }
            // std::cerr << std::endl;

            // for (const auto &data : perf_data) {
            //     double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
            //     std::cerr << std::setw(15) << std::fixed << std::setprecision(2) << mean;
            // }
            // std::cerr << std::endl << std::endl;

        }

        // void perf_init() {
        //     // init perf_header
        //     PerfEvent perf_event;
        //     perf_event.startCounters();
        //     // if (residual_compress_type != "fastpfor" && residual_compress_type != "none")
        //     //     throw std::runtime_error("residual_compress_type not recognised");
        //
        //     perf_event.stopCounters();
        //     std::stringstream header_out;
        //     std::stringstream data_out;
        //     PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
        //     perf_event.printReport(header_out, data_out, 1);
        //     perf_header = split_str(header_out.str(), ',');
        //     perf_data.resize(perf_header.size()); // used for store the performance data
        // }

        // void perf_add(PerfEvent &perf_event) {
        //     std::stringstream header_out;
        //     std::stringstream data_out;
        //     PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
        //     perf_event.printReport(header_out, data_out, 1);
        //     std::vector<double> perf_data_tmp = parseDoubles(data_out.str(), ',');
        //     for (int i = 0; i < perf_data_tmp.size(); i++) {
        //         perf_data[i].push_back(perf_data_tmp[i]);
        //     }
        // }

        void test_model(const std::string input_basename, std::string decode_type) {
            // PerfEvent perf_event;
            // perf_event.startCounters();
            // if (input_basename.back() != '/') {
            //     std::cerr << "Error: output_basename must end with '/'" << std::endl;
            //     return;
            // }
            //
            // perf_event.stopCounters();
            // std::stringstream header_out;
            // std::stringstream data_out;
            // PerfEvent::printCounter(header_out, data_out, "time sec", perf_event.getDuration());
            // perf_event.printReport(header_out, data_out, 1);
            // perf_header = split_str(header_out.str(), ',');
            // perf_data.resize(perf_header.size());

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

            // perf_init();

            uint64_t data_size_tmp = 0;
            in_header.read(reinterpret_cast<char*>(&data_size_tmp), sizeof(data_size_tmp));
            K index_num = 0;
            in_header.read(reinterpret_cast<char*>(&index_num), sizeof(index_num));

            if (epsilon == 0) {
                for (K i = 0; i < index_num; ++i) {
                    size_t partition_size;
                    in_header.read(reinterpret_cast<char*>(&partition_size), sizeof(size_t));

                    std::vector<LICOVariant> partition_index;
                    partition_index.reserve(partition_size);

                    for (size_t j = 0; j < partition_size; ++j) {
                        std::string filename = input_basename + std::to_string(i) + "_" + std::to_string(j) + ".idx";
                        std::ifstream in(filename, std::ios::binary);
                        if (!in) {
                            std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                            throw std::runtime_error("File open error");
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

                        partition_index.push_back(std::move(variant_index));
                        in.close();
                    }

                    lico_enumerator <K> enumerator_tmp = create_enumerator_from_indexes <K>(partition_index);

                    decode_test(enumerator_tmp, decode_type);
                }
            }
            else {
                for (K i = 0; i < index_num; ++i) {
                    std::string filename = input_basename + std::to_string(i) + ".idx";
                    std::ifstream in(filename, std::ios::binary);

                    uint32_t Epsilon_Data;

                    while (in.peek() != EOF) {
                        LICOVariant variant_index;
                        in.read(reinterpret_cast<char*>(&Epsilon_Data), sizeof(Epsilon_Data));
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
                            default: std::cerr << "Unsupported Epsilon Value: " << Epsilon_Data << std::endl; break;
                        }

                        std::visit([&in, Epsilon_Data](auto& index) {
                            index.Epsilon_Data = Epsilon_Data;
                            read_index_data(in, index);
                        }, variant_index);

                        lico_enumerator <K> enumerator_tmp = create_enumerator_from_single_index <K> (variant_index);

                        decode_test(enumerator_tmp, decode_type);
                    }

                    in.close();
                }

            }

            in_header.close();
            result_statistic(decode_type);
        }

    };
}
