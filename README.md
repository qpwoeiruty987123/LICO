# LICO
An SIMD-Aware High-Performance Learned Inverted Index
Compression Framework

This repository provides end-to-end components for building, storing, decoding, and querying LICO-compressed posting lists, with both scalar and AVX-512 SIMD implementation.

## Dataset preparation
The expected dataset directory layout:
```
<data_dir>/
	└─ <DATASET_NAME>/
			 ├─ <DATASET_NAME>.docs        # posting lists in binary
			 ├─ <DATASET_NAME>-2.queries   # one query per line, list IDs separated by spaces
			 ├─ <DATASET_NAME>-3.queries
			 ├─ <DATASET_NAME>-4.queries
			 └─ <DATASET_NAME>-5.queries
```

Binary `.docs` format (little-endian `uint32_t`):
- `data[0] = 1` (magic number, asserted by the builder)
- `data[1] = universe_size`
- Then repeated per-list blocks: `n`, followed by `n` strictly increasing docIDs
	- Each list must be strictly increasing; the max must be < `universe_size` (checked during build)

Query files (`.queries`) format:
- Plain text; each line is a whitespace-separated list of integer list IDs; duplicates on the same line are removed during loading, here is an example of one query:
    
    `2563 876 1236 869 235`
- Use `<DATASET_NAME>-<query_num>.queries` to select 2/3/4/5-term query sets

- Notably, the file `<DATASET_NAME>-5.queries` contains queries that the number of terms is equal of larger than five

The builder runs an internal `data_test` that decodes and checks equality against inputs; any mismatch aborts with details.

## Quickly Start
**Step 1:** You need a server that supports the following instruction sets for compilation: LZCNT, BMI2, and AVX-512 (including AVX512F, AVX512VL, AVX512CD, and AVX512DQ).

**Step 2:** If you want to enable huge pages, comment out line 17 and uncomment line 18 in the CMakeLists.txt file.

**Step 3:** Two convenience bash scripts are provided at the repository root:
- `bdq_lico_greedy.sh`: residuals not FastPFor-compressed (LICO), including build, decode and query test.
- `bdq_lico++_greedy.sh`: residuals FastPFor-compressed (LICO++), including build, decode and query test

**Step 4:** Parameters you should adjust according to your environment in these scripts:
- `result_dir`: The root directory for saving indexes and logs.
- `data_dir`: The dataset root directory containing `<DATASET_NAME>`.
- `code_dir`: The CMake build output directory.
- `dataset`: The name(s) of the dataset(s) to run.
- `read_only`: Set to `t` or `f`. When set to `t`, the script will read an existing index from `index_save_dir` without rebuilding it.


**Step 5:** Run Experiments
- Test LICO
```bash
sh bdq_lico_greedy.sh
```
- Test LICO++
```bash
sh bdq_lico++_greedy.sh
```


Other parameters that use default settings:
- `lico_m`: The upper bound on the number of partitions `m` (limits greedy splitting).
- `compress_type`: Compression method for residuals; options are `none` or `fastpfor`. Note that if you want to test LICO, you must set `compress_type=none` and run the executable with the "none" suffix, and vice versa.
- `decode_type`: `simd` (requires AVX-512) or `normal`
- `query_type`: `AND` (intersection) or `OR` (union)
- `query_num`: 2/3/4/5 (selects the corresponding queries file)


## Architecture overview

Top-level executables (built via CMake from root-level .cpp files):
- `lico_build.cpp`: builds LICO indexes from a dataset and writes them to an output directory.
- `lico_decode.cpp`: loads saved indexes and measures decoding performance (SIMD).
- `lico_query.cpp`: loads selected lists, decodes them, then evaluates list intersection/union performance.

Key headers in `include/` and their roles:
- `config.hpp`: global switches and defaults
    - `HUGEPAGE`: hint to use huge pages when available
	- `RESIDUAL_COMPRESS`: compile-time defines for residual compression
	- `residual_compress_type`: residual codec (currently `fastpfor`)
- `lico_index.hpp`: the core codes to build LICO and LICO++
	- Stores per-segment parametersas bit-packed integer vectors
	- Residuals are stored as bit-packed with magnitude and sign bits (LICO) or ZigZag+FastPFor compressed (LICO++)
- `lico_index_enumerate.hpp`: LICO's high-performance query processing engine
	- `lico_enumerator<K>` holds segments and residuals; implements `normal_decode` and `simd_decode_512i` (AVX-512)
	- Supports `NextGeq` for filter within lists; uses memory re-layout and alignment for SIMD throughput
- `lico_index_build.hpp`: dataset-to-index builder
	- `lico_builder<K, epsilon>` reads `<basename>.docs`, trains/partitions lists, and selects epsilon
		- `epsilon == 0`: partition the list and choose epsilon* per block
		- otherwise: fixed epsilon
	- Includes `statistic_index`, `data_test`, `save_model`, `load_model`
- `lico_index_decode.hpp`: offline decoder/benchmark driver
	- `lico_decoder<K, epsilon>` loads the stored indexes and runs `simd` or `normal` decoding
- `lico_index_query.hpp`: query driver (Intersection/Union)
	- Loads a subset of lists and computes intersection/union using SIMD set ops
- `lico_partition.hpp`: page-based partitioning and epsilon* estimation
	- Greedy/DP partitioning by gap variance; computes error-scaling epsilon* and estimated bit cost per block; merges adjacent blocks with the same epsilon
- `lico_fastpfor.hpp`: FastPFor wrapper for residuals
	- Residuals are ZigZag-transformed and encoded/decoded via `simdfastpfor256`
- `piecewise_linear_model.hpp`: PGM-index-based optimal segmentation, returning LICO’s segment parameters `(delta_y, delta_x, y_b, x_b)`
- `tools.hpp`: utilities
	- Data loading helpers, string parsing, ZigZag encode/decode, and index (de)serialization matching `save_model/load_model`
- `include/sdsl-lite/`: bit-vector utilities used for packing/unpacking
- `external/FastPFor`: the library of [FastPFor](https://github.com/fast-pack/FastPFOR)
- `external/libdivide`: the library of [libdivide](https://github.com/ridiculousfish/libdivide)
- `external/mm_file`: the library of [mm_file](https://github.com/jermp/mm_file)
- `external/Vp2SetOperation`: we reuse the SIMD-based list union algorithm from the [SimSIMD](https://github.com/ashvardanian/simsimd) library, and the SIMD-based list intersection algorithms from the [SIMDsetOperation](https://github.com/tetzank/SIMDSetOperations) and [v2pintersect](https://github.com/mozonaut/vp2intersect) libraries.

Note: your server must support AVX-512 for `SIMD-based` sequential scan, list intersection and list union.

## Building

Standard CMake build from the repository root will produce executables corresponding to the files above. Two variants of each tool are typically produced:
- `*_none`: residuals are not FastPFor-compressed (LICO)
- without suffix: residuals use FastPFor compression (LICO++)


## Notes
- Huge pages: you can define `HUGEPAGE` at compile time to enable huge-page allocation paths (see `CMakeLists.txt`).  You need to ensure that your system has at least **2048 huge pages of 2M size** if you intend to use huge pages.
- Input lists must be strictly increasing and below `universe_size`; this is validated during build
