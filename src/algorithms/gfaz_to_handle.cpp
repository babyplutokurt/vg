#include "gfaz_to_handle.hpp"

#include "../path.hpp"

#include <gfa_compression/decompression_workflow.hpp>
#include <gfa_compression/serialization.hpp>

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <queue>
#include <tuple>

namespace vg {
namespace algorithms {

// Deserialize + decompress a GFAZ file to its in-memory GfaGraph representation.
static GfaGraph load_gfaz_graph(const string& filename, int num_threads) {
    if (filename == "-") {
        throw invalid_argument("GFAZ input from stdin (-) is not supported");
    }

    CompressedData compressed = deserialize_compressed_data(filename);
    GfaGraph graph;
    decompress_gfa(compressed, graph, num_threads);
    return graph;
}

static void write_gfaz_translation(const GFAIDMapInfo& id_map_info, const string& translation_filename) {
    // don't write anything unless we have both an output file and at least one non-trivial mapping
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

static inline size_t segment_count(const GfaGraph& gfaz_graph) {
    return gfaz_graph.node_sequences.empty() ? 0 : gfaz_graph.node_sequences.size() - 1;
}

static inline bool valid_segment_id(const GfaGraph& gfaz_graph, nid_t id) {
    return id > 0 && static_cast<size_t>(id) < gfaz_graph.node_sequences.size();
}

static inline handle_t get_oriented_handle(MutableHandleGraph* graph,
                                           const GfaGraph& gfaz_graph,
                                           NodeId oriented_id) {
    if (oriented_id == 0) {
        throw GFAFormatError("Encountered invalid node id 0 in GFAZ path/walk");
    }
    bool is_reverse = oriented_id < 0;
    nid_t id = is_reverse ? -static_cast<nid_t>(oriented_id) : static_cast<nid_t>(oriented_id);
    if (!valid_segment_id(gfaz_graph, id)) {
        throw GFAFormatError("Encountered out-of-range node id " + to_string(id) + " in GFAZ path/walk");
    }
    return graph->get_handle(id, is_reverse);
}

static void set_numeric_translation(const GfaGraph& gfaz_graph, GFAIDMapInfo* translation) {
    if (!translation) {
        return;
    }
    translation->numeric_mode = true;
    translation->max_id = segment_count(gfaz_graph);
    translation->name_to_id->clear();
    translation->id_to_name.reset();
}

static void add_graph_elements(const GfaGraph& gfaz_graph, MutableHandleGraph* graph) {
    // Nodes are 1-based in GFAZ decompression output.
    for (nid_t id = 1; static_cast<size_t>(id) < gfaz_graph.node_sequences.size(); ++id) {
        graph->create_handle(gfaz_graph.node_sequences[id], id);
    }

    static const string not_blunt = ("error:[gfaz_to_handle_graph] Can only load blunt-ended GFAs. "
        "Try \"bluntifying\" your graph with a tool like <https://github.com/vgteam/GetBlunted>, or "
        "transitively merge overlaps with a pipeline of <https://github.com/ekg/gimbricate> and "
        "<https://github.com/ekg/seqwish>.");

    const auto& links = gfaz_graph.links;
    size_t edge_count = links.from_ids.size();
    for (size_t i = 0; i < edge_count; ++i) {
        nid_t from = links.from_ids[i];
        nid_t to = links.to_ids[i];
        if (!valid_segment_id(gfaz_graph, from) || !valid_segment_id(gfaz_graph, to)) {
            throw GFAFormatError("Encountered out-of-range edge endpoint in GFAZ links");
        }

        uint32_t overlap_num = (i < links.overlap_nums.size()) ? links.overlap_nums[i] : 0;
        char overlap_op = (i < links.overlap_ops.size()) ? links.overlap_ops[i] : '\0';
        bool overlap_is_allowed = (overlap_op == '\0') || (overlap_num == 0 && overlap_op == 'M');
        if (!overlap_is_allowed) {
            string overlap_text = to_string(overlap_num) + overlap_op;
            throw GFAFormatError(not_blunt + " Found edge with a non-null alignment '" + overlap_text + "'.");
        }

        bool from_is_reverse = (i < links.from_orients.size()) ? links.from_orients[i] == '-' : false;
        bool to_is_reverse = (i < links.to_orients.size()) ? links.to_orients[i] == '-' : false;
        graph->create_edge(graph->get_handle(from, from_is_reverse),
                           graph->get_handle(to, to_is_reverse));
    }
}

static unordered_set<string> parse_reference_samples(const GfaGraph& gfaz_graph) {
    unordered_set<string> reference_samples;
    if (gfaz_graph.header_line.empty() || gfaz_graph.header_line[0] != 'H') {
        return reference_samples;
    }

    // Tags are tab-delimited fields on the header line.
    size_t start = 0;
    while (start < gfaz_graph.header_line.size()) {
        size_t end = gfaz_graph.header_line.find('\t', start);
        if (end == string::npos) {
            end = gfaz_graph.header_line.size();
        }
        string tag = gfaz_graph.header_line.substr(start, end - start);
        if (tag.size() >= 5 &&
            std::equal(gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.begin(),
                       gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.end(),
                       tag.begin()) &&
            tag[2] == ':' && tag[3] == 'Z' && tag[4] == ':') {
            reference_samples = gbwtgraph::parse_reference_samples_tag(tag.substr(5));
        }
        start = end + 1;
    }
    return reference_samples;
}

static void add_p_line_paths(const GfaGraph& gfaz_graph,
                             MutablePathMutableHandleGraph* graph,
                             unordered_set<PathSense>* ignore_sense,
                             const unordered_set<string>& reference_samples) {
    size_t count = min(gfaz_graph.path_names.size(), gfaz_graph.paths.size());
    for (size_t i = 0; i < count; ++i) {
        const string& name = gfaz_graph.path_names[i];
        const vector<NodeId>& visits = gfaz_graph.paths[i];

        PathSense sense;
        string sample;
        string locus;
        size_t haplotype;
        size_t phase_block;
        subrange_t subrange;
        PathMetadata::parse_path_name(name,
                                      sense,
                                      sample,
                                      locus,
                                      haplotype,
                                      phase_block,
                                      subrange);

        if (sense == PathSense::HAPLOTYPE && reference_samples.count(sample)) {
            sense = PathSense::REFERENCE;
        } else if (sense == PathSense::REFERENCE &&
                   haplotype != PathMetadata::NO_HAPLOTYPE &&
                   !reference_samples.count(sample)) {
            sense = PathSense::HAPLOTYPE;
            if (phase_block == PathMetadata::NO_PHASE_BLOCK) {
                phase_block = 0;
            }
        }

        if (ignore_sense && ignore_sense->count(sense)) {
            continue;
        }

        string implied_path_name = PathMetadata::create_path_name(sense,
                                                                  sample,
                                                                  locus,
                                                                  haplotype,
                                                                  phase_block,
                                                                  subrange);
        if (graph->has_path(implied_path_name)) {
            throw GFADuplicatePathError(implied_path_name);
        }

        path_handle_t path = graph->create_path(sense,
                                                sample,
                                                locus,
                                                haplotype,
                                                phase_block,
                                                subrange);
        for (NodeId visit : visits) {
            graph->append_step(path, get_oriented_handle(graph, gfaz_graph, visit));
        }
    }
}

static void add_w_line_paths(const GfaGraph& gfaz_graph,
                             MutablePathMutableHandleGraph* graph,
                             const unordered_set<string>& reference_samples) {
    const auto& walks = gfaz_graph.walks;
    for (size_t i = 0; i < walks.walks.size(); ++i) {
        const string sample_name = (i < walks.sample_ids.size()) ? walks.sample_ids[i] : "*";
        const size_t haplotype = (i < walks.hap_indices.size()) ? walks.hap_indices[i] : 0;
        const string contig_name = (i < walks.seq_ids.size()) ? walks.seq_ids[i] : PathMetadata::NO_LOCUS_NAME;
        const int64_t start = (i < walks.seq_starts.size()) ? walks.seq_starts[i] : PathMetadata::NO_END_POSITION;
        const int64_t end = (i < walks.seq_ends.size()) ? walks.seq_ends[i] : PathMetadata::NO_END_POSITION;

        PathSense sense;
        size_t phase_block;
        size_t assigned_haplotype = haplotype;
        string assigned_sample_name;

        if (sample_name == "*") {
            sense = PathSense::GENERIC;
            assigned_sample_name = PathMetadata::NO_SAMPLE_NAME;
            if (assigned_haplotype != 0) {
                throw GFAFormatError("Generic path on omitted (*) sample has nonzero haplotype");
            }
            assigned_haplotype = PathMetadata::NO_HAPLOTYPE;
            phase_block = PathMetadata::NO_PHASE_BLOCK;
        } else {
            if (reference_samples.count(sample_name)) {
                sense = PathSense::REFERENCE;
                phase_block = PathMetadata::NO_PHASE_BLOCK;
            } else {
                sense = PathSense::HAPLOTYPE;
                phase_block = 0;
            }
            assigned_sample_name = sample_name;
        }

        subrange_t assigned_subrange = (start == 0)
            ? PathMetadata::NO_SUBRANGE
            : subrange_t(start, end);

        string implied_path_name = PathMetadata::create_path_name(sense,
                                                                  assigned_sample_name,
                                                                  contig_name,
                                                                  assigned_haplotype,
                                                                  phase_block,
                                                                  assigned_subrange);
        if (graph->has_path(implied_path_name)) {
            throw GFADuplicatePathError(implied_path_name);
        }

        path_handle_t path = graph->create_path(sense,
                                                assigned_sample_name,
                                                contig_name,
                                                assigned_haplotype,
                                                phase_block,
                                                assigned_subrange);
        for (NodeId visit : walks.walks[i]) {
            graph->append_step(path, get_oriented_handle(graph, gfaz_graph, visit));
        }
    }
}

static void add_rgfa_paths(const GfaGraph& gfaz_graph,
                           MutablePathMutableHandleGraph* graph,
                           int64_t max_rgfa_rank) {
    if (max_rgfa_rank < 0) {
        return;
    }

    const OptionalFieldColumn* sn_col = nullptr;
    const OptionalFieldColumn* so_col = nullptr;
    const OptionalFieldColumn* sr_col = nullptr;
    for (const auto& col : gfaz_graph.segment_optional_fields) {
        if (col.tag == "SN" && col.type == 'Z') {
            sn_col = &col;
        } else if (col.tag == "SO" && col.type == 'i') {
            so_col = &col;
        } else if (col.tag == "SR" && col.type == 'i') {
            sr_col = &col;
        }
    }
    if (!sn_col || !so_col || !sr_col) {
        return;
    }

    vector<size_t> sn_offsets(sn_col->string_lengths.size() + 1, 0);
    for (size_t i = 0; i < sn_col->string_lengths.size(); ++i) {
        sn_offsets[i + 1] = sn_offsets[i] + sn_col->string_lengths[i];
    }

    struct RGFAVisit {
        int64_t offset;
        nid_t node_id;
        size_t length;
    };
    struct RGFAPath {
        int64_t rank = 0;
        bool rank_set = false;
        vector<RGFAVisit> visits;
    };
    unordered_map<string, RGFAPath> by_path;

    const int64_t missing_i64 = numeric_limits<int64_t>::min();
    for (nid_t node_id = 1; static_cast<size_t>(node_id) < gfaz_graph.node_sequences.size(); ++node_id) {
        size_t idx = node_id - 1;
        if (idx >= sn_col->string_lengths.size() ||
            idx >= so_col->int_values.size() ||
            idx >= sr_col->int_values.size()) {
            continue;
        }
        uint32_t sn_len = sn_col->string_lengths[idx];
        int64_t so = so_col->int_values[idx];
        int64_t sr = sr_col->int_values[idx];
        if (sn_len == 0 || so == missing_i64 || sr == missing_i64 || sr > max_rgfa_rank) {
            continue;
        }

        string path_name = sn_col->concatenated_strings.substr(sn_offsets[idx], sn_len);
        auto& path_info = by_path[path_name];
        if (path_info.rank_set && path_info.rank != sr) {
            throw GFAFormatError("rGFA path " + path_name + " has conflicting ranks " +
                                 to_string(sr) + " and " + to_string(path_info.rank));
        }
        path_info.rank = sr;
        path_info.rank_set = true;
        path_info.visits.push_back({so, node_id, gfaz_graph.node_sequences[node_id].size()});
    }

    for (auto& kv : by_path) {
        const string& path_name = kv.first;
        auto& visits = kv.second.visits;
        sort(visits.begin(), visits.end(), [](const RGFAVisit& a, const RGFAVisit& b) {
            return a.offset < b.offset;
        });

        path_handle_t current_path;
        int64_t next_expected = -1;
        bool have_path = false;

        for (const auto& visit : visits) {
            if (!have_path || visit.offset != next_expected) {
                subrange_t subrange = (visit.offset == 0)
                    ? PathMetadata::NO_SUBRANGE
                    : subrange_t(visit.offset, PathMetadata::NO_END_POSITION);
                string implied_path_name = PathMetadata::create_path_name(PathSense::GENERIC,
                                                                          PathMetadata::NO_SAMPLE_NAME,
                                                                          path_name,
                                                                          PathMetadata::NO_HAPLOTYPE,
                                                                          PathMetadata::NO_PHASE_BLOCK,
                                                                          subrange);
                if (graph->has_path(implied_path_name)) {
                    throw GFADuplicatePathError(implied_path_name);
                }
                current_path = graph->create_path(PathSense::GENERIC,
                                                  PathMetadata::NO_SAMPLE_NAME,
                                                  path_name,
                                                  PathMetadata::NO_HAPLOTYPE,
                                                  PathMetadata::NO_PHASE_BLOCK,
                                                  subrange);
                have_path = true;
            }
            graph->append_step(current_path, graph->get_handle(visit.node_id, false));
            next_expected = visit.offset + visit.length;
        }
    }
}

void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          GFAIDMapInfo* translation,
                          int num_threads) {
    GfaGraph gfaz_graph = load_gfaz_graph(filename, num_threads);
    set_numeric_translation(gfaz_graph, translation);
    add_graph_elements(gfaz_graph, graph);
}

void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          const string& translation_filename,
                          int num_threads) {
    GFAIDMapInfo id_map_info;
    gfaz_to_handle_graph(filename, graph, &id_map_info, num_threads);
    write_gfaz_translation(id_map_info, translation_filename);
}

void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               GFAIDMapInfo* translation,
                               int64_t max_rgfa_rank,
                               unordered_set<PathSense>* ignore_sense,
                               int num_threads) {
    GfaGraph gfaz_graph = load_gfaz_graph(filename, num_threads);
    set_numeric_translation(gfaz_graph, translation);
    add_graph_elements(gfaz_graph, graph);

    unordered_set<string> reference_samples = parse_reference_samples(gfaz_graph);
    add_p_line_paths(gfaz_graph, graph, ignore_sense, reference_samples);
    add_w_line_paths(gfaz_graph, graph, reference_samples);
    add_rgfa_paths(gfaz_graph, graph, max_rgfa_rank);
}

void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               int64_t max_rgfa_rank,
                               const string& translation_filename,
                               unordered_set<PathSense>* ignore_sense,
                               int num_threads) {
    GFAIDMapInfo id_map_info;
    gfaz_to_path_handle_graph(filename, graph, &id_map_info, max_rgfa_rank, ignore_sense, num_threads);
    write_gfaz_translation(id_map_info, translation_filename);
}

}
}
