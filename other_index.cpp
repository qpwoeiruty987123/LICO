#include <other_index.hpp>

#include "lico_index_build.hpp"

std::string log_filename = "";
std::string function_type = "";
std::string query_filename = "";
std::string query_type = "";

void build_collection_others(const std::string &input_basename, const std::string &index_type, const std::string &output_filename) {
    typedef other_sequence::other_index_build OTHER_INDEX;
    OTHER_INDEX other_index;
    other_index.build_index(input_basename + ".docs", index_type);
    // other_index.load_model(output_filename);
    other_index.statistic_index(log_filename);
    // other_index.decode_test();
    other_index.data_test(input_basename + ".docs", index_type);
    other_index.save_model(output_filename);
}

void decode_collection_others(const std::string &input_filename) {
    typedef other_sequence::other_index_decode OTHER_INDEX;
    OTHER_INDEX other_index;
    other_index.load_model(input_filename);
    other_index.result_statistics();
}

void query_collection_others(const std::string &input_basename, const std::string &output_filename) {
    typedef other_sequence::other_index_query OTHER_INDEX;
    OTHER_INDEX index;
    index.test_query(output_filename, query_filename, query_type, input_basename + ".docs");
}

template <uint64_t epsilon=64>
void build_collection_la_vector(const std::string input_basename, const std::string output_basename) {
    typedef lico_sequence::lico_builder<uint32_t, epsilon> PGM_INDEX_BUILDER;
    PGM_INDEX_BUILDER index;
    index.build_model(input_basename + ".docs", 1);
    index.data_test(input_basename + ".docs");
    index.save_model(output_basename);
    if(!log_filename.empty()) { // save covered and residual
        index.statistic_index(log_filename, true);
    }
}


int main(int argc, const char** argv)
{

    int mandatory = 7;
    if (argc < mandatory) {
        std::cerr << "Usage: " << argv[0] << ":\n" << "\t <index_type> <collection_basename> <index_filename> <function_type> <query_filename> <query_type> <log_filename> <epsilon>" << std::endl;
        return 1;
    }

    const std::string index_type = argv[1];
    const std::string input_basename = argv[2];
    const std::string output_basename = argv[3];
    function_type = argv[4];
    query_filename = argv[5];
    query_type = argv[6];
    log_filename = argv[7];
    const std::string epsilonstr = argv[8];
    const uint64_t epsilon = static_cast<uint64_t>(std::stoi(epsilonstr));

    if (index_type=="fastpfor") {
        if (function_type == "build")
            build_collection_others(input_basename, index_type, output_basename);
        else if (function_type == "decode")
            decode_collection_others(output_basename);
        else if (function_type == "query")
            query_collection_others(input_basename, output_basename);
    } else if (index_type=="la_vector") {
        if (function_type == "build") {
            switch (epsilon) {
                case 15: build_collection_la_vector<15>(input_basename, output_basename); break;
                case 31: build_collection_la_vector<31>(input_basename, output_basename); break;
                case 63: build_collection_la_vector<63>(input_basename, output_basename); break;
                case 127: build_collection_la_vector<127>(input_basename, output_basename); break;
                case 255: build_collection_la_vector<255>(input_basename, output_basename); break;
                case 511: build_collection_la_vector<511>(input_basename, output_basename); break;
                case 1023: build_collection_la_vector<1023>(input_basename, output_basename); break;
                case 2047: build_collection_la_vector<2047>(input_basename, output_basename); break;
                case 4095: build_collection_la_vector<4095>(input_basename, output_basename); break;
                case 8191: build_collection_la_vector<8191>(input_basename, output_basename); break;
                case 16383: build_collection_la_vector<16383>(input_basename, output_basename); break;
                case 32767: build_collection_la_vector<32767>(input_basename, output_basename); break;
                case 65535: build_collection_la_vector<65535>(input_basename, output_basename); break;
                case 131071: build_collection_la_vector<131071>(input_basename, output_basename); break;
                case 262143: build_collection_la_vector<262143>(input_basename, output_basename); break;
                default: std::cerr << "Unsupported Epsilon Value: " << epsilon << std::endl; break;
            }
        }
        else if (function_type == "decode") {

        }
            // decode_collection_others(output_basename);
        else if (function_type == "query") {

        }
            // query_collection_others(input_basename, output_basename);
    }



    return 0;
}
