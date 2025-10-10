#pragma once
#include <string>

# ifndef  RESIDUAL_COMPRESS
# define RESIDUAL_COMPRESS 0
# endif

static constexpr size_t block_size_rice = 8192;

static std::string residual_compress_type = "fastpfor";

static bool use_huge_pages = true;

