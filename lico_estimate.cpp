#include <lico_epsilon_estimation.hpp>

void build_collection_lico(const std::string input_basename, const std::string collection_name) {
    typedef lico_sequence::lico_estimater<uint32_t> PGM_INDEX_ESTIMATER;
    PGM_INDEX_ESTIMATER index;
    if (collection_name == "uniform-1M")
        index.test_epsilon_uniform();
    else
        index.test_epsilon(input_basename + ".docs");
}

int main(int argc, const char** argv)
{
    std::ios::sync_with_stdio(0);
    int mandatory = 2;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <collection_path> <collection_name>" << std::endl;
        return 1;
    }

    const std::string input_basename = argv[1];
    const std::string collection_name = argv[2];

    build_collection_lico(input_basename, collection_name);

    return 0;
}
