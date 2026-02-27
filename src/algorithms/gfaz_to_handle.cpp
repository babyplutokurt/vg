#include "gfaz_to_handle.hpp"

#include "../utility.hpp"

#include <gfa_compression/decompression_workflow.hpp>
#include <gfa_compression/gfa_writer.hpp>
#include <gfa_compression/serialization.hpp>

#include <fstream>

namespace vg {
namespace algorithms {

// Materialize a decompressed temporary GFA and return its filename.
static string materialize_gfaz_as_gfa(const string& filename, int num_threads) {
    if (filename == "-") {
        throw invalid_argument("GFAZ input from stdin (-) is not supported");
    }

    CompressedData compressed = deserialize_compressed_data(filename);
    GfaGraph decompressed;
    decompress_gfa(compressed, decompressed, num_threads);

    string temp_gfa = temp_file::create("vg-gfaz-input-");
    write_gfa(decompressed, temp_gfa);
    return temp_gfa;
}

void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          GFAIDMapInfo* translation,
                          int num_threads) {
    string temp_gfa = materialize_gfaz_as_gfa(filename, num_threads);
    try {
        gfa_to_handle_graph(temp_gfa, graph, translation);
        temp_file::remove(temp_gfa);
    } catch (...) {
        temp_file::remove(temp_gfa);
        throw;
    }
}

void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          const string& translation_filename,
                          int num_threads) {
    GFAIDMapInfo id_map_info;
    gfaz_to_handle_graph(filename, graph, &id_map_info, num_threads);
    // Reuse GFA import translation writer behavior by going through the existing overload.
    if (!translation_filename.empty() && !id_map_info.numeric_mode) {
        ofstream trans_file(translation_filename);
        if (!trans_file) {
            throw runtime_error("error:[gfaz_to_handle_graph] Unable to open output translation file: " + translation_filename);
        }
        for (const auto& mapping : *id_map_info.name_to_id) {
            trans_file << "T\t" << mapping.first << "\t" << mapping.second << "\n";
        }
    }
}

void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               GFAIDMapInfo* translation,
                               int64_t max_rgfa_rank,
                               unordered_set<PathSense>* ignore_sense,
                               int num_threads) {
    string temp_gfa = materialize_gfaz_as_gfa(filename, num_threads);
    try {
        gfa_to_path_handle_graph(temp_gfa, graph, translation, max_rgfa_rank, ignore_sense);
        temp_file::remove(temp_gfa);
    } catch (...) {
        temp_file::remove(temp_gfa);
        throw;
    }
}

void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               int64_t max_rgfa_rank,
                               const string& translation_filename,
                               unordered_set<PathSense>* ignore_sense,
                               int num_threads) {
    GFAIDMapInfo id_map_info;
    gfaz_to_path_handle_graph(filename, graph, &id_map_info, max_rgfa_rank, ignore_sense, num_threads);
    if (!translation_filename.empty() && !id_map_info.numeric_mode) {
        ofstream trans_file(translation_filename);
        if (!trans_file) {
            throw runtime_error("error:[gfaz_to_path_handle_graph] Unable to open output translation file: " + translation_filename);
        }
        for (const auto& mapping : *id_map_info.name_to_id) {
            trans_file << "T\t" << mapping.first << "\t" << mapping.second << "\n";
        }
    }
}

}
}
