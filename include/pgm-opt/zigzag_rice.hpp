# pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>

// ---------- sign <-> zigzag ----------
/**
 * 将有符号整数x转换为zigzag编码（无符号）。
 */
inline uint32_t zigzag(int32_t x) {
    return (static_cast<uint32_t>(x) << 1) ^ static_cast<uint32_t>(x >> 31);
}

/**
 * 将zigzag编码的无符号整数u还原为有符号整数。
 */
inline int32_t unzigzag(uint32_t u) {
    return static_cast<int32_t>((u >> 1) ^ -(static_cast<int32_t>(u & 1)));
}

// ---------- tiny bitstream ----------
/**
 * BitWriter：用于逐位写入数据并输出为字节流。
 */
class BitWriter {
public:
    BitWriter() : buf(0), nbits(0) {}

    /**
     * 写入n位的整数v到缓冲区。
     */
    void write_bits(uint32_t v, int n) {
        buf = (buf << n) | (v & ((1U << n) - 1));
        nbits += n;
        while (nbits >= 8) {
            nbits -= 8;
            uint8_t byte = (buf >> nbits) & 0xFF;
            out.push_back(byte);
        }
    }

    /**
     * 写入unary编码：q个1后跟一个0
     */
    void write_unary(uint32_t q) {
        while (q >= 32) {
            write_bits((1U << 32) - 1, 32);
            q -= 32;
        }
        if (q)
            write_bits((1U << q) - 1, q);
        write_bits(0, 1);
    }

    /**
     * 输出所有缓冲内容为字节流，最后不足8位补零
     */
    std::vector<uint8_t> finish() {
        if (nbits)
            write_bits(0, 8 - nbits); // pad zeros
        return out;
    }

private:
    uint64_t buf;           // 用于缓冲写入的比特
    int nbits;              // 缓冲区内的比特数量
    std::vector<uint8_t> out; // 输出字节数组
};

/**
 * BitReader：用于逐位读取字节流数据。
 */
class BitReader {
public:
    BitReader(const uint8_t* data, size_t size)
        : data(data), data_size(size), index(0), buf(0), nbits(0) {}

    /**
     * 读取n位比特为整数。
     */
    uint32_t read_bits(int n) {
        while (nbits < n) {
            if (index >= data_size)
                throw std::runtime_error("bitstream underflow");
            buf = (buf << 8) | data[index++];
            nbits += 8;
        }
        nbits -= n;
        return (buf >> nbits) & ((1U << n) - 1);
    }

    /**
     * 读取unary编码，返回1的数量
     */
    uint32_t read_unary() {
        uint32_t q = 0;
        while (true) {
            uint32_t b = read_bits(1);
            if (b == 0)
                return q;
            q += 1;
        }
    }

private:
    const uint8_t* data;    // 输入数据指针
    size_t data_size;       // 输入数据长度
    size_t index;           // 当前读取索引
    uint64_t buf;           // 缓冲区
    int nbits;              // 缓冲区内比特
};

// ---------- Rice core ----------
/**
 * 估算Rice编码块所需的总比特数
 */
size_t rice_estimated_bits(const std::vector<uint32_t>& block_u, int k) {
    size_t total = 0;
    for (size_t i = 0; i < block_u.size(); ++i)
        total += (block_u[i] >> k) + 1 + k;
    return total;
}

/**
 * 自动选择最佳k（0~kmax），使编码比特最少
 */
int choose_k(const std::vector<uint32_t>& block_u, int kmax = 12) {
    int best_k = 0;
    size_t best_cost = rice_estimated_bits(block_u, 0);
    for (int k = 1; k <= kmax; ++k) {
        size_t cost = rice_estimated_bits(block_u, k);
        if (cost < best_cost) {
            best_k = k;
            best_cost = cost;
        }
    }
    return best_k;
}

/**
 * 编码一个区块到BitWriter
 */
void rice_encode_block(BitWriter& bw, const std::vector<uint32_t>& block_u, int k) {
    // 用4位写入k参数
    bw.write_bits(k, 4);
    for (size_t i = 0; i < block_u.size(); ++i) {
        uint32_t u = block_u[i];
        uint32_t q = u >> k;
        uint32_t r = u & ((1U << k) - 1);
        bw.write_unary(q);
        if (k)
            bw.write_bits(r, k);
    }
}

/**
 * 解码一个区块并返回结果
 */
std::vector<uint32_t> rice_decode_block(BitReader& br, size_t n) {
    int k = br.read_bits(4);
    std::vector<uint32_t> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t q = br.read_unary();
        uint32_t r = k ? br.read_bits(k) : 0;
        out.push_back((q << k) | r);
    }
    return out;
}

static size_t all_block_size = 8192;
// ---------- Public API ----------
/**
 * 压缩残差整数，返回字节流
 */
std::vector<uint8_t> compress_residuals(const std::vector<int32_t>& ints, size_t block_size = all_block_size) {
    // 1) zigzag编码
    std::vector<uint32_t> U;
    U.reserve(ints.size());
    for (size_t i = 0; i < ints.size(); ++i)
        U.emplace_back(zigzag(ints[i]));

    BitWriter bw;
    // 写入总数量（32位）
    bw.write_bits(static_cast<uint32_t>(U.size()), 32);
    // 块式Rice编码
    for (size_t i = 0; i < U.size(); i += block_size) {
        size_t end = std::min(U.size(), i + block_size);
        std::vector<uint32_t> block(U.begin() + i, U.begin() + end);
        int k = choose_k(block);
        rice_encode_block(bw, block, k);
    }
    return bw.finish();
}

/**
 * 解压残差整数，返回原始有符号整数
 */
std::vector<int32_t> decompress_residuals(const uint8_t* data, size_t data_size, size_t block_size = all_block_size) {
    BitReader br(data, data_size);
    uint32_t n = br.read_bits(32);
    std::vector<uint32_t> out_u;
    out_u.reserve(n);
    size_t remaining = n;
    while (remaining > 0) {
        size_t take = std::min(block_size, remaining);
        std::vector<uint32_t> block = rice_decode_block(br, take);
        out_u.insert(out_u.end(), block.begin(), block.end());
        remaining -= take;
    }
    // Unzigzag
    std::vector<int32_t> out;
    out.reserve(n);
    for (size_t i = 0; i < out_u.size(); ++i)
        out.emplace_back(unzigzag(out_u[i]));
    return out;
}