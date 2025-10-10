# pragma once
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <../config.hpp>

/**
 * BitWriter：write bits and output bytes
 */
class BitWriter {
public:
    BitWriter() : buf(0), nbits(0) {}


    void write_bits(uint32_t v, int n) {
        buf = (buf << n) | (v & ((1U << n) - 1));
        nbits += n;
        while (nbits >= 8) {
            nbits -= 8;
            uint8_t byte = (buf >> nbits) & 0xFF;
            out.push_back(byte);
        }
    }


    void write_unary(uint32_t q) {
        while (q >= 32) {
            write_bits((1U << 32) - 1, 32);
            q -= 32;
        }
        if (q)
            write_bits((1U << q) - 1, q);
        write_bits(0, 1);
    }


    std::vector<uint8_t> finish() {
        if (nbits)
            write_bits(0, 8 - nbits); // pad zeros
        return out;
    }

private:
    uint64_t buf;
    int nbits;
    std::vector<uint8_t> out;
};

/**
 * BitReader：read bits
 */
class BitReader {
public:
    BitReader() : buf(0), nbits(0) {}
    BitReader(const uint8_t* data, size_t size)
        : data(data), data_size(size), index(0), buf(0), nbits(0) {}

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
    const uint8_t* data;
    size_t data_size;
    size_t index;
    uint64_t buf;
    int nbits;
};

// ---------- Rice core ----------
size_t rice_estimated_bits(const std::vector<uint32_t>& block_u, int k) {
    size_t total = 0;
    for (size_t i = 0; i < block_u.size(); ++i)
        total += (block_u[i] >> k) + 1 + k;
    return total;
}

/**
 * automatically choose k to minimize the total cost
 */
int choose_k(const std::vector<uint32_t>& block_u, int kmax = 12) { // 可以修改
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

/* rice block encode */
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

/* rice block decode */
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


// ---------- Public API ----------
std::vector<uint8_t> compress_residuals_rice(const std::vector<int32_t>& ints) {
    // 1) zigzag
    std::vector<uint32_t> U;
    U.reserve(ints.size());
    for (size_t i = 0; i < ints.size(); ++i)
        U.push_back(zigzag(ints[i]));

    BitWriter bw;
    bw.write_bits(static_cast<uint32_t>(U.size()), 32);
    // rice block encode
    for (size_t i = 0; i < U.size(); i += block_size_rice) {
        size_t end = std::min(U.size(), i + block_size_rice);
        std::vector<uint32_t> block(U.begin() + i, U.begin() + end);
        int k = choose_k(block);
        rice_encode_block(bw, block, k);
    }
    return bw.finish();
}

std::vector<int32_t> decompress_residuals_rice(const uint8_t* data, size_t data_size) {
    BitReader br(data, data_size);
    uint32_t n = br.read_bits(32);
    std::vector<uint32_t> out_u;
    out_u.reserve(n);
    size_t remaining = n;
    while (remaining > 0) {
        size_t take = std::min(block_size_rice, remaining);
        std::vector<uint32_t> block = rice_decode_block(br, take);
        out_u.insert(out_u.end(), block.begin(), block.end());
        remaining -= take;
    }
    // Unzigzag
    std::vector<int32_t> out;
    out.reserve(n);

    for (size_t i = 0; i < out_u.size(); ++i)
        out.push_back(unzigzag(out_u[i]));

    return out;
}


class RiceDecoder {
public:
    RiceDecoder() = default;

    RiceDecoder(const uint8_t* data, size_t data_size)
        : remaining(0)
        , block_size(8192) {
        br = new BitReader(data, data_size);
        remaining = br->read_bits(32);
    }

    uint32_t decompress_residuals_block(int* out) {
        uint32_t take = std::min(block_size, remaining);
        if (take == 0) return 0;
        int k = br -> read_bits(4);
        for (uint32_t i = 0; i < take; ++i) {
            uint32_t q = br -> read_unary();
            uint32_t r = k ? br -> read_bits(k) : 0;
            out[i] = (q << k) | r;
            out[i] = unzigzag(out[i]);
        }
        remaining -= take;
        return take;
    }

    uint32_t block_size;
    BitReader* br;
    uint32_t remaining;
};