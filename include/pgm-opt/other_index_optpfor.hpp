#pragma once

#include <fstream>
#include <../../external/mm_file/include/mm_file/mm_file.hpp>
#include <iostream>
#include <vector>
#include <chrono>

// #include "codecfactory.h"
// #include "optpfor.h"
// #include "simple16.h"
#include  <optpfor.h>


typedef uint32_t K;

namespace other_sequence {

    #define PREFIX_SUM                               \
    output[0] += cur_base;                          \
    for (uint32_t k = 1; k != block_size; ++k) { \
    output[k] += output[k - 1] + 1;                \
    }

    // 前向声明
    class other_index;

    static size_t block_size = 128;
    static constexpr uint32_t block_size_u32 = 128;

    // 单个列表的 FastPFor 编码器/解码器
    class fastpfor_list {
    public:
        // 构造函数：可选预分配缓冲区大小
        fastpfor_list() : compressed_data_(), n(0) {}

        // 编码：输入原始数据，压缩并存储
        void encode(std::vector<K>& data) {
            if (data.empty()) {
                compressed_data_.clear();
                block_header_.clear();
                compress_size_ = 0;
                n = 0;
                return;
            }

            n = data.size();
            compressed_data_.resize(n + 4096);

            K* compressed_data_begin = compressed_data_.data();
            K* compressed_data_pointer = compressed_data_.data();

            const uint64_t block_num = (data.size() + block_size - 1) / block_size;

            block_header_.resize(block_num * 2 - 1); // n block_max, n-1 block_endpoint

            // generate gap vector
            std::vector<K> docs_buf(block_size * 2);
            K* docs_it = data.data();
            code_type codec_;
            auto start = std::chrono::high_resolution_clock::now();
            uint32_t last_doc(-1);
            for (uint64_t b = 0; b < block_num; b++) {
                size_t cur_block_size = std::min<size_t>(block_size, n - b * block_size);
                size_t true_compress_size = docs_buf.size();
                for (size_t i = 0; i < cur_block_size; ++i) {
                    K doc(*docs_it++);
                    docs_buf[i] = doc - last_doc - 1;
                    last_doc = doc;
                }
                block_header_[b] = last_doc; // max_doc of each block

                if (cur_block_size == block_size)
                    codec_.encodeBlock(docs_buf.data(), compressed_data_pointer, true_compress_size);
                // else
                    // codec_.encodeArray(docs_buf.data(), cur_block_size, compressed_data_pointer, true_compress_size);
                else {
                    // std::cerr << "Block Size: " << cur_block_size << std::endl;
                    // codec_.encodeArray(docs_buf.data(), cur_block_size, compressed_data_pointer, true_compress_size);
                    for (size_t i = 0; i < cur_block_size; ++i) {
                        compressed_data_pointer[i] = docs_buf[i];
                    }
                    true_compress_size = cur_block_size;
                }
                compressed_data_pointer += true_compress_size;

                if (b != block_num - 1) {
                    block_header_[block_num + b] = compressed_data_pointer - compressed_data_begin; // endpoint
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            total_duration_ = duration.count();

            compressed_data_.resize(compressed_data_pointer - compressed_data_begin);
            compressed_data_.shrink_to_fit();
            compress_size_ = compressed_data_.size() + block_header_.size();
        }

        void decode(K* output) {
            const uint64_t block_num = (n + block_size - 1) / block_size;
            // assert(block_num > 0);

            uint32_t* compressed_data_begin = compressed_data_.data();
            uint32_t* block_max_begin = block_header_.data();
            uint32_t* block_endpoint_begin = block_max_begin + block_num;


            code_type codec_;
            uint32_t endpoint = 0;
            uint32_t i = 0;


            auto start = std::chrono::high_resolution_clock::now();
            for (; i != block_num - 1; ++i) {
                uint32_t cur_base = (i ? block_max_begin[i - 1] : uint32_t(-1)) + 1;
                uint32_t const* ptr = compressed_data_begin + endpoint;

                // codec_.decodeArray(ptr, block_endpoint_begin[i] - endpoint, output, block_size);
                codec_.decodeBlock(ptr, output, block_size);

                output[0] += cur_base;
                for (int k = 1; k < 128; ++k) {
                    output[k] += output[k - 1] + 1;
                }

                // PREFIX_SUM

                endpoint = block_endpoint_begin[i];
                output += block_size;
            }

            // last block
            uint32_t cur_base = (i ? block_max_begin[i - 1] : uint32_t(-1)) + 1;
            uint32_t const* ptr = compressed_data_begin + endpoint;
            uint32_t last_block_size = n - (block_num - 1) * block_size;
            if (last_block_size == block_size)
                codec_.decodeBlock(ptr, output, block_size);
            else {
                for (size_t i = 0; i < last_block_size; ++i) {
                    output[i] = ptr[i];
                }
            }
                // codec_.decodeArray(ptr, compressed_data_.size() - endpoint, output, last_block_size);

            output[0] += cur_base;
            for (uint32_t k = 1; k < last_block_size; ++k) {
                output[k] += output[k - 1] + 1;
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            total_duration_ = duration.count();
        }

        // 获取压缩后数据（用于保存）
        const std::vector<uint32_t>& get_compressed_data() const {
            return compressed_data_;
        }

        const std::vector<uint32_t>& get_header_data() const {
            return block_header_;
        }

        // 从外部数据恢复压缩内容（用于加载）
        void set_compressed_data(std::vector<uint32_t> &data, size_t original_size) {
            compressed_data_ = std::move(data);
            n = original_size;
        }

        void set_header_data(std::vector<uint32_t> &data) {
            block_header_ = std::move(data);
            compress_size_ = compressed_data_.size() + block_header_.size(); // 同步
        }

        size_t size() const {
            return (compressed_data_.size() + block_header_.size()) * sizeof(uint32_t);
        }


        size_t original_size() const {
            return n;
        }

        bool empty() const {
            return compressed_data_.empty();
        }

        size_t time_cost_ns() const {
            return total_duration_;
        }

    private:
        std::vector<uint32_t> compressed_data_;
        std::vector<uint32_t> block_header_;
        size_t n;
        size_t compress_size_;
        size_t total_duration_ = 0;

        // 共享 codec 实例（静态或每个实例持有一个）
        typedef FastPForLib::OPTPFor<4, FastPForLib::Simple16<false>> code_type ;
        // FastPForLib::CODECFactory factory_;
        // FastPForLib::IntegerCODEC& codec_ = *factory_.getFromName("optpfor");
    };


    // builder
    class other_index_build {
    private:
        uint64_t data_size = 0;
        uint64_t data_compressed_size = 0;

        // 存储每个 list 的压缩对象
        std::vector<fastpfor_list> encoded_lists;
        std::vector<uint64_t> list_lengths;  // 原始长度，用于解码分配空间

        bool use_u8 = false;  // 可扩展支持其他编码（如 leco）

    public:
        void build_index(const std::string& input_basename, const std::string& index_type) {
            std::cerr << "Read File [Build]: " << input_basename << std::endl;
            mm::file_source<uint32_t> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;

            K list_count = 0;
            for (size_t i = 2; i < input.size();) {
                uint32_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);
                data_size += n;
                list_lengths.push_back(n);
                list_count++;
                if (list_count % 100 == 0)
                    std::cerr << list_count << " ";
                if (index_type == "fastpfor" || index_type == "optpfor") {
                    fastpfor_list list_encoder;
                    list_encoder.encode(sequence);
                    encoded_lists.push_back(std::move(list_encoder));
                }
                // else if (index_type == "leco") { ... }

                i += n + 1;
            }
        }

        void statistic_index(std::string output_basename = "") {
            if (output_basename.empty()) {
                std::cerr << "Warning: No output basename provided. Using stdout for statistics." << std::endl;
            }

            std::ofstream file;
            if (!output_basename.empty()) {
                file.open(output_basename + ".statistic_log.txt");
                if (!file.is_open()) {
                    std::cerr << "Error: Cannot open " << output_basename << ".statistic_log.txt for writing." << std::endl;
                    return;
                }
                // 输出文件头
                file << "#list_id,original_elements,original_bytes,compressed_bytes,compression_ratio(%)\n";
            }

            // 全局统计
            uint64_t total_original_bytes = 0;  // 原始总字节数
            uint64_t total_compressed_bytes = 0; // 压缩后总字节数
            uint64_t total_duration = 0;

            // 遍历每个列表进行统计
            for (size_t i = 0; i < encoded_lists.size(); ++i) {
                const auto& list_encoder = encoded_lists[i];
                const size_t original_size = list_encoder.original_size(); // 元素个数
                const size_t compressed_bytes = list_encoder.size();       // 压缩后字节数

                const size_t original_bytes = original_size * sizeof(K);   // 原始字节数

                total_original_bytes += original_bytes;
                total_compressed_bytes += compressed_bytes;
                total_duration += list_encoder.time_cost_ns();

                double compression_ratio = 0.0;
                if (original_bytes > 0) {
                    compression_ratio = static_cast<double>(original_bytes) / compressed_bytes * 100.0;
                }

                // 写入每个 list 的详细信息
                if (file.is_open()) {
                    file << i << " "
                         << original_size << " "
                         << original_bytes << " "
                         << compressed_bytes << " "
                         << compression_ratio << "%\n";
                }
            }

            // 全局指标计算
            double global_compression_ratio = 0.0;
            double bits_per_integer = 0.0;

            if (total_original_bytes > 0) {
                global_compression_ratio = static_cast<double>(total_original_bytes) /  total_compressed_bytes * 100.0;
            }

            if (data_size > 0) {  // data_size 是原始整数总个数
                bits_per_integer = static_cast<double>(total_compressed_bytes * 8) / data_size;
            }

            // 输出全局统计
            std::ostream& out = file.is_open() ? static_cast<std::ostream&>(file) : std::cout;

            out << "\n# Global Statistics\n";
            out << "Total original size (elements): " << data_size << "\n";
            out << "Total original bytes: " << total_original_bytes << " Bytes\n";
            out << "Total compressed bytes: " << total_compressed_bytes << " Bytes, " << static_cast<long double> (total_compressed_bytes) / 1024 / 1024 / 1024 <<  " GiB \n";
            out << "Global compression ratio: " << global_compression_ratio << "%\n";
            out << "Compression factor (original/compressed): " << (total_compressed_bytes > 0 ? static_cast<double>(total_original_bytes) / total_compressed_bytes : 0.0) << "%\n";
            out << "Average bits per integer (BPI): " << bits_per_integer << " bits\n";

            std::cerr << "\nBuild Time: " << static_cast<long double> (total_duration) / 1e9 << "\n";
            std::cerr << "\nTotal compressed bytes: " << total_compressed_bytes << " Bytes, " << static_cast<long double> (total_compressed_bytes) / 1024 / 1024 / 1024 << " GiB, Global compression ratio: " << global_compression_ratio << "%, " << "Average bits per integer (BPI): " << bits_per_integer << " bits\n";

            if (file.is_open()) {
                file.close();
                std::cerr << "Statistics saved to: " << output_basename << ".statistic_log.txt" << std::endl;
            }
        }

        void data_test(const std::string& input_basename, const std::string& index_type) {
            std::cerr << "Read File [Test]: " << input_basename << std::endl;
            mm::file_source<uint32_t> input(input_basename.c_str(), mm::advice::sequential);
            K const* data = input.data();
            std::cerr << "Universe Size: " << data[1] << std::endl;

            size_t data_unequal_test = 0;
            K idx = 0;
            for (size_t i = 2; i < input.size();) {
                uint32_t n = data[i];
                std::vector<K> sequence(data + i + 1, data + i + n + 1);

                std::vector<K> result(sequence.size());

                encoded_lists[idx].decode(result.data());

                assert(result.size() == sequence.size());

                for (auto j = 0; j < result.size(); j++) {
                    if (sequence[j] != result[j]) {
                        data_unequal_test++;
                        // std::cerr << "Unequal Value: " << result_decode[j] << " " << sequence[j] << " " << posi << " " << j << std::endl;
                    }
                }

                i += n + 1;
                idx++;
            }
            std::cerr << "Unequal postings: " << data_unequal_test << std::endl << std::endl;
        }

        void decode_test() {
            uint64_t total_duration = 0;
            uint64_t total_integers = 0;
            for (uint32_t i = 0; i < encoded_lists.size(); i++) {
                auto &list = encoded_lists[i];
                total_integers += list.original_size();
                std::vector<K> result(list.original_size());
                for (int32_t j = 0; j < 3; j++)
                    for (int32_t k = 0; k < list.original_size(); k++)
                        result[k] = j + k;

                auto start = std::chrono::high_resolution_clock::now();
                list.decode(result.data());
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                total_duration += duration.count();
            }
            std::cerr << "Total duration: " << total_duration << std::endl;
            std::cerr << "Decode per int: " << static_cast<long double> (total_duration) / total_integers << " ns" << std::endl;
        }

        uint64_t index_size() {
            data_compressed_size = 0;
            for (const auto& lst : encoded_lists) {
                data_compressed_size += lst.size();
            }
            return data_compressed_size;
        }

        void save_model(const std::string& output_basename) {
            if (output_basename.empty() || output_basename.back() != '/') {
                std::cerr << "Error: output_basename must end with '/'" << std::endl;
                return;
            }

            std::cerr << "Save index to: " << output_basename << std::endl;

            std::ofstream out_header(output_basename + "idx.size", std::ios::binary);
            if (!out_header) {
                std::cerr << "Error: Cannot open idx.size for writing." << std::endl;
                return;
            }

            out_header.write(reinterpret_cast<const char*>(&data_size), sizeof(data_size));
            use_u8 = false; // 当前只支持 fastpfor (uint32)
            out_header.write(reinterpret_cast<const char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size = encoded_lists.size();
            out_header.write(reinterpret_cast<const char*>(&total_list_size), sizeof(size_t));

            for (size_t i = 0; i < total_list_size; ++i) {
                const auto& compressed_data = encoded_lists[i].get_compressed_data();
                const auto& block_header = encoded_lists[i].get_header_data();
                size_t comp_size = compressed_data.size();
                size_t header_size = block_header.size();
                size_t orig_size = encoded_lists[i].original_size();

                std::string filename = output_basename + std::to_string(i) + ".idx";
                std::ofstream out(filename, std::ios::binary);
                if (!out) {
                    std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
                    throw std::runtime_error("file open error");
                }

                out.write(reinterpret_cast<const char*>(&orig_size), sizeof(orig_size));
                out.write(reinterpret_cast<const char*>(&comp_size), sizeof(comp_size));
                out.write(reinterpret_cast<const char*>(&header_size), sizeof(header_size));
                out.write(reinterpret_cast<const char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                out.write(reinterpret_cast<const char*>(block_header.data()), header_size * sizeof(uint32_t));
                out.close();
            }
            out_header.close();
        }

        void load_model(const std::string& input_basename) {
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
            in_header.read(reinterpret_cast<char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size;
            in_header.read(reinterpret_cast<char*>(&total_list_size), sizeof(size_t));
            in_header.close();

            encoded_lists.clear();
            encoded_lists.reserve(total_list_size);
            list_lengths.clear();
            list_lengths.reserve(total_list_size);

            for (size_t i = 0; i < total_list_size; ++i) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    continue;
                }

                size_t orig_size, comp_size, header_size;
                in.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
                in.read(reinterpret_cast<char*>(&comp_size), sizeof(comp_size));
                in.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

                std::vector<uint32_t> compressed_data(comp_size);
                in.read(reinterpret_cast<char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                std::vector<uint32_t> block_header(header_size);
                in.read(reinterpret_cast<char*>(block_header.data()), header_size * sizeof(uint32_t));
                in.close();

                fastpfor_list list;
                list.set_compressed_data(compressed_data, orig_size);
                list.set_header_data(block_header);
                encoded_lists.push_back(std::move(list));
                list_lengths.push_back(orig_size);
            }
        }

        // 示例：获取第 i 个列表的解码结果
        // std::vector<K> decode_list(size_t idx) {
        //     if (idx >= encoded_lists.size()) {
        //         return {};
        //     }
        //     return encoded_lists[idx].decode();
        // }

        // 示例：解码到已有数组
        void decode_list(size_t idx, K* output) {
            if (idx < encoded_lists.size()) {
                encoded_lists[idx].decode(output);
            }
        }
    };

    // decode
    class other_index_decode {
    private:
        uint64_t data_size = 0;
        uint64_t total_duration = 0;
        uint64_t total_integers = 0;
        // 存储每个 list 的压缩对象

        std::vector<uint64_t> list_lengths;  // 原始长度，用于解码分配空间

        bool use_u8 = false;  // 可扩展支持其他编码（如 leco）

    public:
        void decode_test(fastpfor_list &list) {
            total_integers += list.original_size() + 5;
            std::vector<K> result(list.original_size() );

            for (int32_t j = 0; j < 3; j++)
                for (int32_t k = 0; k < list.original_size(); k++)
                    result[k] = j + k;

            // auto start = std::chrono::high_resolution_clock::now();
            list.decode(result.data());
            // auto end = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            // total_duration += duration.count();
            total_duration += list.time_cost_ns();
        }

        void result_statistics() {
            std::cerr << "Total duration: " << total_duration << std::endl;
            std::cerr << "Decode per int: " << static_cast<long double> (total_duration) / total_integers << " ns" << std::endl;
        }

        void load_model(const std::string& input_basename) {
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
            in_header.read(reinterpret_cast<char*>(&use_u8), sizeof(use_u8));

            size_t total_list_size;
            in_header.read(reinterpret_cast<char*>(&total_list_size), sizeof(size_t));
            in_header.close();


            for (size_t i = 0; i < total_list_size; ++i) {
                std::string filename = input_basename + std::to_string(i) + ".idx";
                std::ifstream in(filename, std::ios::binary);
                if (!in) {
                    std::cerr << "Error: Cannot open " << filename << " for reading." << std::endl;
                    continue;
                }

                size_t orig_size, comp_size, header_size;
                in.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
                in.read(reinterpret_cast<char*>(&comp_size), sizeof(comp_size));
                in.read(reinterpret_cast<char*>(&header_size), sizeof(header_size));

                std::vector<uint32_t> compressed_data(comp_size);
                in.read(reinterpret_cast<char*>(compressed_data.data()), comp_size * sizeof(uint32_t));
                std::vector<uint32_t> block_header(header_size);
                in.read(reinterpret_cast<char*>(block_header.data()), header_size * sizeof(uint32_t));
                in.close();

                fastpfor_list list;
                list.set_compressed_data(compressed_data, orig_size);
                list.set_header_data(block_header);
                decode_test(list);
            }
        }
    };

} // namespace other_sequence