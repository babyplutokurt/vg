#include "gfa_to_handle.hpp"

#include "../path.hpp"

#include <GFAz/codec/codec.hpp>
#include <GFAz/codec/serialization.hpp>
#include <GFAz/io/gfa_write_utils.hpp>

#include <gbwtgraph/utils.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <tuple>

namespace vg {
namespace algorithms {
using gfaz::gfa_write_utils::SequenceOffsets;
using gfaz::gfa_write_utils::build_offsets;
using gfaz::gfa_write_utils::decode_rules;
using gfaz::gfa_write_utils::decompress_optional_column;
using gfaz::gfa_write_utils::decompress_string_column;
using gfaz::CompressedData;
using gfaz::LinkData;
using gfaz::NodeId;
using gfaz::OptionalFieldColumn;
using gfaz::GFAZ_MAGIC;
using gfaz::deserialize_compressed_data;
namespace Codec = gfaz::Codec;

struct StreamingGFAZPaths {
  string header_line;
  size_t node_count = 0;
  vector<size_t> node_lengths;
  vector<OptionalFieldColumn> segment_optional_fields;

  LinkData links;

  vector<int32_t> rules_first;
  vector<int32_t> rules_second;
  vector<int32_t> paths_flat;
  vector<int32_t> walks_flat;
  SequenceOffsets path_offsets;
  SequenceOffsets walk_offsets;
  SequenceOffsets original_path_offsets;
  SequenceOffsets original_walk_offsets;

  vector<string> path_names;
  vector<string> walk_sample_ids;
  vector<uint32_t> walk_hap_indices;
  vector<string> walk_seq_ids;
  vector<int64_t> walk_seq_starts;
  vector<int64_t> walk_seq_ends;

  uint32_t min_rule_id = 0;
};

bool GFAzParser::looks_like_gfaz(const string& filename) {
  if (filename.empty() || filename == "-") {
    return false;
  }
  ifstream in(filename, ios::binary);
  if (!in) {
    return false;
  }
  uint32_t magic = 0;
  in.read(reinterpret_cast<char*>(&magic), sizeof(uint32_t));
  if (in.gcount() < static_cast<streamsize>(sizeof(uint32_t))) {
    return false;
  }
  return magic == GFAZ_MAGIC;
}

static CompressedData load_gfaz_compressed(const string &filename) {
  if (filename == "-") {
    throw invalid_argument("GFAZ input from stdin (-) is not supported");
  }
  return deserialize_compressed_data(filename);
}

static void write_gfaz_translation(const GFAIDMapInfo &id_map_info,
                                   const string &translation_filename) {
  if (!translation_filename.empty() && !id_map_info.numeric_mode) {
    ofstream trans_file(translation_filename);
    if (!trans_file) {
      throw runtime_error("error:[gfaz_to_handle_graph] Unable to open output "
                          "translation file: " +
                          translation_filename);
    }
    for (const auto &mapping : *id_map_info.name_to_id) {
      trans_file << "T\t" << mapping.first << "\t" << mapping.second << "\n";
    }
  }
}

static inline bool valid_segment_id(size_t segment_count, nid_t id) {
  return id > 0 && static_cast<size_t>(id) <= segment_count;
}

static inline handle_t get_oriented_handle(MutableHandleGraph *graph,
                                           size_t segment_count,
                                           NodeId oriented_id) {
  if (oriented_id == 0) {
    throw GFAFormatError("Encountered invalid node id 0 in GFAZ path/walk");
  }
  bool is_reverse = oriented_id < 0;
  nid_t id = is_reverse ? -static_cast<nid_t>(oriented_id)
                        : static_cast<nid_t>(oriented_id);
  if (!valid_segment_id(segment_count, id)) {
    throw GFAFormatError("Encountered out-of-range node id " + to_string(id) +
                         " in GFAZ path/walk");
  }
  return graph->get_handle(id, is_reverse);
}

static void set_numeric_translation(size_t segment_count,
                                    GFAIDMapInfo *translation) {
  if (!translation) {
    return;
  }
  translation->numeric_mode = true;
  translation->max_id = segment_count;
  translation->name_to_id->clear();
  translation->id_to_name.reset();
}

static void add_edge_elements(const LinkData &links, size_t node_count,
                              MutableHandleGraph *graph) {
  static const string not_blunt =
      ("error:[gfaz_to_handle_graph] Can only load blunt-ended GFAs. "
       "Try \"bluntifying\" your graph with a tool like "
       "<https://github.com/vgteam/GetBlunted>, or "
       "transitively merge overlaps with a pipeline of "
       "<https://github.com/ekg/gimbricate> and "
       "<https://github.com/ekg/seqwish>.");

  size_t edge_count = links.from_ids.size();
  for (size_t i = 0; i < edge_count; ++i) {
    nid_t from = links.from_ids[i];
    nid_t to = links.to_ids[i];
    if (!valid_segment_id(node_count, from) || !valid_segment_id(node_count, to)) {
      throw GFAFormatError(
          "Encountered out-of-range edge endpoint in GFAZ links");
    }

    uint32_t overlap_num =
        (i < links.overlap_nums.size()) ? links.overlap_nums[i] : 0;
    char overlap_op =
        (i < links.overlap_ops.size()) ? links.overlap_ops[i] : '\0';
    bool overlap_is_allowed =
        (overlap_op == '\0') || (overlap_num == 0 && overlap_op == 'M');
    if (!overlap_is_allowed) {
      string overlap_text = to_string(overlap_num) + overlap_op;
      throw GFAFormatError(not_blunt +
                           " Found edge with a non-null alignment '" +
                           overlap_text + "'.");
    }

    bool from_is_reverse =
        (i < links.from_orients.size()) ? links.from_orients[i] == '-' : false;
    bool to_is_reverse =
        (i < links.to_orients.size()) ? links.to_orients[i] == '-' : false;
    graph->create_edge(graph->get_handle(from, from_is_reverse),
                       graph->get_handle(to, to_is_reverse));
  }
}

static unordered_set<string> parse_reference_samples(const string &header_line) {
  unordered_set<string> reference_samples;
  if (header_line.empty() || header_line[0] != 'H') {
    return reference_samples;
  }

  size_t start = 0;
  while (start < header_line.size()) {
    size_t end = header_line.find('\t', start);
    if (end == string::npos) {
      end = header_line.size();
    }
    string tag = header_line.substr(start, end - start);
    if (tag.size() >= 5 &&
        std::equal(gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.begin(),
                   gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG.end(),
                   tag.begin()) &&
        tag[2] == ':' && tag[3] == 'Z' && tag[4] == ':') {
      auto parsed_reference_samples =
          gbwtgraph::parse_reference_samples_tag(tag.substr(5));
      reference_samples.clear();
      reference_samples.insert(parsed_reference_samples.begin(),
                               parsed_reference_samples.end());
    }
    start = end + 1;
  }
  return reference_samples;
}

static void expand_rule(uint32_t rule_id, bool reverse,
                        const vector<int32_t> &first,
                        const vector<int32_t> &second, uint32_t min_id,
                        uint32_t max_id, vector<NodeId> &out) {
  const uint32_t idx = rule_id - min_id;
  const int32_t a = first[idx];
  const int32_t b = second[idx];

  if (!reverse) {
    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_id && abs_a < max_id) {
      expand_rule(abs_a, a < 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(a);
    }

    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_id && abs_b < max_id) {
      expand_rule(abs_b, b < 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(b);
    }
  } else {
    const uint32_t abs_b = static_cast<uint32_t>(std::abs(b));
    if (abs_b >= min_id && abs_b < max_id) {
      expand_rule(abs_b, b >= 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(-b);
    }

    const uint32_t abs_a = static_cast<uint32_t>(std::abs(a));
    if (abs_a >= min_id && abs_a < max_id) {
      expand_rule(abs_a, a >= 0, first, second, min_id, max_id, out);
    } else {
      out.push_back(-a);
    }
  }
}

static vector<NodeId> decode_sequence_at_index(
    const vector<int32_t> &flat, const SequenceOffsets &compressed_offsets,
    const SequenceOffsets &original_offsets, size_t index,
    const vector<int32_t> &rules_first, const vector<int32_t> &rules_second,
    uint32_t min_rule_id, int delta_round) {
  if (index + 1 >= compressed_offsets.size()) {
    throw out_of_range("GFAZ traversal index out of range");
  }

  const size_t start = compressed_offsets[index];
  const size_t end = compressed_offsets[index + 1];
  if (end > flat.size()) {
    throw runtime_error("GFAZ traversal block is truncated");
  }

  const uint32_t max_rule_id =
      min_rule_id + static_cast<uint32_t>(rules_first.size());
  const size_t original_length =
      (index + 1 < original_offsets.size())
          ? (original_offsets[index + 1] - original_offsets[index])
          : (end - start);

  vector<NodeId> decoded;
  decoded.reserve(original_length);
  for (size_t pos = start; pos < end; ++pos) {
    const NodeId node = flat[pos];
    const uint32_t abs_id = static_cast<uint32_t>(std::abs(node));
    if (abs_id >= min_rule_id && abs_id < max_rule_id) {
      expand_rule(abs_id, node < 0, rules_first, rules_second, min_rule_id,
                  max_rule_id, decoded);
    } else {
      decoded.push_back(node);
    }
  }

  vector<vector<NodeId>> seqs(1);
  seqs[0] = std::move(decoded);
  for (int round = 0; round < delta_round; ++round) {
    Codec::inverse_delta_transform(seqs);
  }
  return std::move(seqs[0]);
}

static void add_p_line_paths_streaming(
    const StreamingGFAZPaths &gfaz_paths, MutablePathMutableHandleGraph *graph,
    unordered_set<PathSense> *ignore_sense,
    const unordered_set<string> &reference_samples, int delta_round) {
  size_t encoded_path_count =
      gfaz_paths.path_offsets.empty() ? 0 : gfaz_paths.path_offsets.size() - 1;
  size_t count = min(gfaz_paths.path_names.size(), encoded_path_count);
  for (size_t i = 0; i < count; ++i) {
    const string &name = gfaz_paths.path_names[i];
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.paths_flat, gfaz_paths.path_offsets,
        gfaz_paths.original_path_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, delta_round);

    PathSense sense;
    string sample;
    string locus;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(name, sense, sample, locus, haplotype,
                                  phase_block, subrange);

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

    string implied_path_name = PathMetadata::create_path_name(
        sense, sample, locus, haplotype, phase_block, subrange);
    if (graph->has_path(implied_path_name)) {
      throw GFADuplicatePathError(implied_path_name);
    }

    path_handle_t path = graph->create_path(sense, sample, locus, haplotype,
                                            phase_block, subrange);
    for (NodeId visit : visits) {
      graph->append_step(path,
                         get_oriented_handle(graph, gfaz_paths.node_count, visit));
    }
  }
}

static void add_w_line_paths_streaming(
    const StreamingGFAZPaths &gfaz_paths, MutablePathMutableHandleGraph *graph,
    const unordered_set<string> &reference_samples, int delta_round) {
  size_t count = gfaz_paths.walk_offsets.size();
  if (count > 0) {
    --count;
  }
  for (size_t i = 0; i < count; ++i) {
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.walks_flat, gfaz_paths.walk_offsets,
        gfaz_paths.original_walk_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, delta_round);
    const string sample_name =
        (i < gfaz_paths.walk_sample_ids.size()) ? gfaz_paths.walk_sample_ids[i]
                                                : "*";
    const size_t haplotype =
        (i < gfaz_paths.walk_hap_indices.size()) ? gfaz_paths.walk_hap_indices[i]
                                                 : 0;
    const string contig_name = (i < gfaz_paths.walk_seq_ids.size())
                                   ? gfaz_paths.walk_seq_ids[i]
                                   : PathMetadata::NO_LOCUS_NAME;
    const int64_t start = (i < gfaz_paths.walk_seq_starts.size())
                              ? gfaz_paths.walk_seq_starts[i]
                              : PathMetadata::NO_END_POSITION;
    const int64_t end = (i < gfaz_paths.walk_seq_ends.size())
                            ? gfaz_paths.walk_seq_ends[i]
                            : PathMetadata::NO_END_POSITION;

    PathSense sense;
    size_t phase_block;
    size_t assigned_haplotype = haplotype;
    string assigned_sample_name;

    if (sample_name == "*") {
      sense = PathSense::GENERIC;
      assigned_sample_name = PathMetadata::NO_SAMPLE_NAME;
      if (assigned_haplotype != 0) {
        throw GFAFormatError(
            "Generic path on omitted (*) sample has nonzero haplotype");
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

    subrange_t assigned_subrange =
        (start == 0) ? PathMetadata::NO_SUBRANGE : subrange_t(start, end);

    string implied_path_name = PathMetadata::create_path_name(
        sense, assigned_sample_name, contig_name, assigned_haplotype,
        phase_block, assigned_subrange);
    if (graph->has_path(implied_path_name)) {
      throw GFADuplicatePathError(implied_path_name);
    }

    path_handle_t path =
        graph->create_path(sense, assigned_sample_name, contig_name,
                           assigned_haplotype, phase_block, assigned_subrange);
    for (NodeId visit : visits) {
      graph->append_step(path,
                         get_oriented_handle(graph, gfaz_paths.node_count, visit));
    }
  }
}

static void dispatch_rgfa_visits(
    size_t node_count, const vector<OptionalFieldColumn>& segment_optional_fields,
    const vector<size_t>& node_lengths, int64_t max_rgfa_rank,
    const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback);

static void add_rgfa_paths(size_t node_count,
                           const vector<OptionalFieldColumn> &segment_optional_fields,
                           MutablePathMutableHandleGraph *graph,
                           const vector<size_t> &node_lengths,
                           int64_t max_rgfa_rank) {
  if (max_rgfa_rank < 0) {
    return;
  }

  using rgfa_cache_t = unordered_map<string, pair<path_handle_t, int64_t>>;
  rgfa_cache_t rgfa_cache;
  dispatch_rgfa_visits(
      node_count, segment_optional_fields, node_lengths, max_rgfa_rank,
      [&](nid_t id, int64_t offset, size_t length, const string &path_name,
          int64_t path_rank) {
        (void)path_rank;
        auto found = rgfa_cache.find(path_name);
        if (found != rgfa_cache.end() && found->second.second != offset) {
          rgfa_cache.erase(found);
          found = rgfa_cache.end();
        }
        if (found == rgfa_cache.end()) {
          subrange_t subrange =
              (offset == 0)
                  ? PathMetadata::NO_SUBRANGE
                  : subrange_t(offset, PathMetadata::NO_END_POSITION);
          string implied_path_name = PathMetadata::create_path_name(
              PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, path_name,
              PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK, subrange);
          if (graph->has_path(implied_path_name)) {
            throw GFADuplicatePathError(implied_path_name);
          }
          path_handle_t path = graph->create_path(
              PathSense::GENERIC, PathMetadata::NO_SAMPLE_NAME, path_name,
              PathMetadata::NO_HAPLOTYPE, PathMetadata::NO_PHASE_BLOCK, subrange);
          found = rgfa_cache.emplace_hint(found, path_name,
                                          std::make_pair(path, offset));
        }

        graph->append_step(found->second.first, graph->get_handle(id, false));
        found->second.second += length;
      });
}


static GFAParser::chars_t as_chars(const string& value) {
  return make_pair(value.begin(), value.end());
}

static GFAParser::tag_list_t collect_optional_tags_for_row(
    const vector<OptionalFieldColumn>& columns,
    vector<size_t>& string_offsets,
    vector<size_t>& byte_offsets,
    size_t row_index) {
  GFAParser::tag_list_t tags;
  tags.reserve(columns.size());
  for (size_t i = 0; i < columns.size(); ++i) {
    const auto& col = columns[i];
    switch (col.type) {
      case 'A':
        if (row_index < col.char_values.size()) {
          tags.push_back(col.tag + ":A:" + string(1, col.char_values[row_index]));
        }
        break;
      case 'i':
        if (row_index < col.int_values.size()) {
          tags.push_back(col.tag + ":i:" + to_string(col.int_values[row_index]));
        }
        break;
      case 'f':
        if (row_index < col.float_values.size()) {
          tags.push_back(col.tag + ":f:" + to_string(col.float_values[row_index]));
        }
        break;
      case 'Z':
      case 'J':
      case 'H':
        if (row_index < col.string_lengths.size()) {
          uint32_t len = col.string_lengths[row_index];
          string value = col.concatenated_strings.substr(string_offsets[i], len);
          string_offsets[i] += len;
          tags.push_back(col.tag + ":" + string(1, col.type) + ":" + value);
        }
        break;
      case 'B':
        if (row_index < col.b_lengths.size() && row_index < col.b_subtypes.size()) {
          uint32_t len = col.b_lengths[row_index];
          size_t offset = byte_offsets[i];
          byte_offsets[i] += len;
          string value = col.tag + ":B:" + string(1, col.b_subtypes[row_index]);
          for (uint32_t j = 0; j < len && offset + j < col.b_concat_bytes.size(); ++j) {
            value.push_back(',');
            value += to_string(col.b_concat_bytes[offset + j]);
          }
          tags.push_back(std::move(value));
        }
        break;
      default:
        break;
    }
  }
  return tags;
}

void gfaz_to_handle_graph(const string &filename, MutableHandleGraph *graph,
                          GFAIDMapInfo *translation) {
  CompressedData compressed = load_gfaz_compressed(filename);

  string segment_sequences =
      Codec::zstd_decompress_string(compressed.segment_sequences_zstd);
  vector<uint32_t> segment_lengths =
      Codec::zstd_decompress_uint32_vector(compressed.segment_seq_lengths_zstd);
  const size_t node_count = segment_lengths.size();

  set_numeric_translation(node_count, translation);

  size_t segment_seq_offset = 0;
  for (size_t i = 0; i < segment_lengths.size(); ++i) {
    nid_t id = static_cast<nid_t>(i + 1);
    uint32_t len = segment_lengths[i];
    if (segment_seq_offset + len > segment_sequences.size()) {
      throw runtime_error("GFAZ segment sequence column is truncated");
    }
    graph->create_handle(segment_sequences.substr(segment_seq_offset, len), id);
    segment_seq_offset += len;
  }

  LinkData links;
  links.from_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_from_ids_zstd, compressed.num_links);
  links.to_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_to_ids_zstd, compressed.num_links);
  links.from_orients = Codec::decompress_orientations(
      compressed.link_from_orients_zstd, compressed.num_links);
  links.to_orients = Codec::decompress_orientations(
      compressed.link_to_orients_zstd, compressed.num_links);
  links.overlap_nums =
      Codec::zstd_decompress_uint32_vector(compressed.link_overlap_nums_zstd);
  links.overlap_ops =
      Codec::zstd_decompress_char_vector(compressed.link_overlap_ops_zstd);
  add_edge_elements(links, node_count, graph);
}

void gfaz_to_handle_graph(const string &filename, MutableHandleGraph *graph,
                          const string &translation_filename) {
  GFAIDMapInfo id_map_info;
  gfaz_to_handle_graph(filename, graph, &id_map_info);
  write_gfaz_translation(id_map_info, translation_filename);
}

void gfaz_to_path_handle_graph(const string &filename,
                               MutablePathMutableHandleGraph *graph,
                               GFAIDMapInfo *translation,
                               int64_t max_rgfa_rank,
                               unordered_set<PathSense> *ignore_sense) {
  CompressedData compressed = load_gfaz_compressed(filename);
  StreamingGFAZPaths gfaz_paths;

  {
    string segment_sequences =
        Codec::zstd_decompress_string(compressed.segment_sequences_zstd);
    vector<uint32_t> segment_lengths =
        Codec::zstd_decompress_uint32_vector(compressed.segment_seq_lengths_zstd);
    gfaz_paths.header_line = compressed.header_line;
    gfaz_paths.node_count = segment_lengths.size();
    gfaz_paths.node_lengths.resize(gfaz_paths.node_count + 1, 0);

    size_t segment_seq_offset = 0;
    for (size_t i = 0; i < segment_lengths.size(); ++i) {
      nid_t id = static_cast<nid_t>(i + 1);
      uint32_t len = segment_lengths[i];
      if (segment_seq_offset + len > segment_sequences.size()) {
        throw runtime_error("GFAZ segment sequence column is truncated");
      }
      graph->create_handle(segment_sequences.substr(segment_seq_offset, len), id);
      gfaz_paths.node_lengths[id] = len;
      segment_seq_offset += len;
    }
  }

  set_numeric_translation(gfaz_paths.node_count, translation);

  gfaz_paths.links.from_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_from_ids_zstd, compressed.num_links);
  gfaz_paths.links.to_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_to_ids_zstd, compressed.num_links);
  gfaz_paths.links.from_orients = Codec::decompress_orientations(
      compressed.link_from_orients_zstd, compressed.num_links);
  gfaz_paths.links.to_orients = Codec::decompress_orientations(
      compressed.link_to_orients_zstd, compressed.num_links);
  gfaz_paths.links.overlap_nums =
      Codec::zstd_decompress_uint32_vector(compressed.link_overlap_nums_zstd);
  gfaz_paths.links.overlap_ops =
      Codec::zstd_decompress_char_vector(compressed.link_overlap_ops_zstd);
  add_edge_elements(gfaz_paths.links, gfaz_paths.node_count, graph);

  if (max_rgfa_rank >= 0) {
    gfaz_paths.segment_optional_fields.reserve(
        compressed.segment_optional_fields_zstd.size());
    for (const auto &column : compressed.segment_optional_fields_zstd) {
      gfaz_paths.segment_optional_fields.push_back(
          decompress_optional_column(column));
    }
  }

  gfaz_paths.path_names =
      decompress_string_column(compressed.names_zstd, compressed.name_lengths_zstd);
  gfaz_paths.walk_sample_ids = decompress_string_column(
      compressed.walk_sample_ids_zstd, compressed.walk_sample_id_lengths_zstd);
  gfaz_paths.walk_hap_indices =
      Codec::zstd_decompress_uint32_vector(compressed.walk_hap_indices_zstd);
  gfaz_paths.walk_seq_ids = decompress_string_column(
      compressed.walk_seq_ids_zstd, compressed.walk_seq_id_lengths_zstd);
  gfaz_paths.walk_seq_starts = Codec::decompress_varint_int64(
      compressed.walk_seq_starts_zstd, compressed.walk_lengths.size());
  gfaz_paths.walk_seq_ends = Codec::decompress_varint_int64(
      compressed.walk_seq_ends_zstd, compressed.walk_lengths.size());

  auto decoded_rules = decode_rules(compressed);
  gfaz_paths.rules_first = std::move(decoded_rules.first);
  gfaz_paths.rules_second = std::move(decoded_rules.second);
  gfaz_paths.min_rule_id = compressed.min_rule_id();

  if (!compressed.paths_zstd.payload.empty()) {
    gfaz_paths.paths_flat = Codec::zstd_decompress_int32_vector(
        compressed.paths_zstd);
  }
  if (!compressed.walks_zstd.payload.empty()) {
    gfaz_paths.walks_flat = Codec::zstd_decompress_int32_vector(
        compressed.walks_zstd);
  }
  gfaz_paths.path_offsets = build_offsets(compressed.sequence_lengths);
  gfaz_paths.walk_offsets = build_offsets(compressed.walk_lengths);
  gfaz_paths.original_path_offsets =
      build_offsets(compressed.original_path_lengths);
  gfaz_paths.original_walk_offsets =
      build_offsets(compressed.original_walk_lengths);

  unordered_set<string> reference_samples =
      parse_reference_samples(gfaz_paths.header_line);
  add_p_line_paths_streaming(gfaz_paths, graph, ignore_sense, reference_samples,
                             compressed.delta_round);
  add_w_line_paths_streaming(gfaz_paths, graph, reference_samples,
                             compressed.delta_round);
  add_rgfa_paths(gfaz_paths.node_count, gfaz_paths.segment_optional_fields, graph,
                 gfaz_paths.node_lengths, max_rgfa_rank);
}

void gfaz_to_path_handle_graph(const string &filename,
                               MutablePathMutableHandleGraph *graph,
                               int64_t max_rgfa_rank,
                               const string &translation_filename,
                               unordered_set<PathSense> *ignore_sense) {
  GFAIDMapInfo id_map_info;
  gfaz_to_path_handle_graph(filename, graph, &id_map_info, max_rgfa_rank,
                            ignore_sense);
  write_gfaz_translation(id_map_info, translation_filename);
}

static GFAParser::visit_iteratee_t make_gfaz_visit_iteratee(const vector<NodeId>& visits) {
  return [visits](const GFAParser::visit_step_t& visit_step) {
    if (visits.empty()) {
      visit_step(-1, string(), false);
      return;
    }
    for (size_t i = 0; i < visits.size(); ++i) {
      NodeId visit = visits[i];
      bool is_reverse = visit < 0;
      nid_t node_id = is_reverse ? -static_cast<nid_t>(visit) : static_cast<nid_t>(visit);
      if (!visit_step(i, to_string(node_id), is_reverse)) {
        return;
      }
    }
  };
}

static void dispatch_rgfa_visits(size_t node_count,
                                 const vector<OptionalFieldColumn>& segment_optional_fields,
                                 const vector<size_t>& node_lengths,
                                 int64_t max_rgfa_rank,
                                 const function<void(nid_t, int64_t, size_t, const string&, int64_t)>& callback) {
  if (max_rgfa_rank < 0) {
    return;
  }

  const OptionalFieldColumn* sn_col = nullptr;
  const OptionalFieldColumn* so_col = nullptr;
  const OptionalFieldColumn* sr_col = nullptr;
  for (const auto& col : segment_optional_fields) {
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

  for (nid_t node_id = 1; static_cast<size_t>(node_id) <= node_count; ++node_id) {
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
      throw GFAFormatError("rGFA path " + path_name +
                           " has conflicting ranks " + to_string(sr) + " and " +
                           to_string(path_info.rank));
    }
    path_info.rank = sr;
    path_info.rank_set = true;
    size_t length = (static_cast<size_t>(node_id) < node_lengths.size()) ? node_lengths[node_id] : 0;
    path_info.visits.push_back({so, node_id, length});
  }

  for (auto& kv : by_path) {
    auto& path_name = kv.first;
    auto& path_info = kv.second;
    sort(path_info.visits.begin(), path_info.visits.end(), [](const RGFAVisit& a, const RGFAVisit& b) {
      return a.offset < b.offset;
    });
    for (const auto& visit : path_info.visits) {
      callback(visit.node_id, visit.offset, visit.length, path_name, path_info.rank);
    }
  }
}

void GFAzParser::parse(istream& in) {
  (void)in;
  throw invalid_argument("GFAZ input from streams is not supported directly");
}

void GFAzParser::parse(const string& filename) {
  CompressedData compressed = load_gfaz_compressed(filename);
  StreamingGFAZPaths gfaz_paths;

  {
    string segment_sequences =
        Codec::zstd_decompress_string(compressed.segment_sequences_zstd);
    vector<uint32_t> segment_lengths =
        Codec::zstd_decompress_uint32_vector(compressed.segment_seq_lengths_zstd);
    gfaz_paths.header_line = compressed.header_line;
    gfaz_paths.node_count = segment_lengths.size();
    gfaz_paths.node_lengths.resize(gfaz_paths.node_count + 1, 0);

    auto& map_info = this->id_map();
    map_info.numeric_mode = true;
    map_info.max_id = gfaz_paths.node_count;
    map_info.name_to_id->clear();
    map_info.id_to_name.reset();
    for (nid_t id = 1; static_cast<size_t>(id) <= gfaz_paths.node_count; ++id) {
      map_info.name_to_id->emplace(to_string(id), id);
    }

    gfaz_paths.segment_optional_fields.reserve(compressed.segment_optional_fields_zstd.size());
    for (const auto& column : compressed.segment_optional_fields_zstd) {
      gfaz_paths.segment_optional_fields.push_back(decompress_optional_column(column));
    }

    if (!gfaz_paths.header_line.empty() && gfaz_paths.header_line[0] == 'H') {
      auto parsed = GFAParser::parse_h(gfaz_paths.header_line);
      for (auto& listener : this->header_listeners) {
        listener(get<0>(parsed));
      }
    }

    vector<size_t> optional_string_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    vector<size_t> optional_byte_offsets(gfaz_paths.segment_optional_fields.size(), 0);
    size_t segment_seq_offset = 0;
    for (size_t i = 0; i < segment_lengths.size(); ++i) {
      nid_t id = static_cast<nid_t>(i + 1);
      uint32_t len = segment_lengths[i];
      if (segment_seq_offset + len > segment_sequences.size()) {
        throw runtime_error("GFAZ segment sequence column is truncated");
      }
      string sequence = segment_sequences.substr(segment_seq_offset, len);
      auto tags = collect_optional_tags_for_row(gfaz_paths.segment_optional_fields,
                                                optional_string_offsets,
                                                optional_byte_offsets,
                                                i);
      for (auto& listener : this->node_listeners) {
        listener(id, as_chars(sequence), tags);
      }
      gfaz_paths.node_lengths[id] = len;
      segment_seq_offset += len;
    }
  }

  gfaz_paths.links.from_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_from_ids_zstd, compressed.num_links);
  gfaz_paths.links.to_ids = Codec::decompress_delta_varint_uint32(
      compressed.link_to_ids_zstd, compressed.num_links);
  gfaz_paths.links.from_orients = Codec::decompress_orientations(
      compressed.link_from_orients_zstd, compressed.num_links);
  gfaz_paths.links.to_orients = Codec::decompress_orientations(
      compressed.link_to_orients_zstd, compressed.num_links);
  gfaz_paths.links.overlap_nums =
      Codec::zstd_decompress_uint32_vector(compressed.link_overlap_nums_zstd);
  gfaz_paths.links.overlap_ops =
      Codec::zstd_decompress_char_vector(compressed.link_overlap_ops_zstd);

  for (size_t i = 0; i < gfaz_paths.links.from_ids.size(); ++i) {
    string overlap;
    if (i < gfaz_paths.links.overlap_ops.size() && gfaz_paths.links.overlap_ops[i] != '\0') {
      overlap = to_string(i < gfaz_paths.links.overlap_nums.size() ? gfaz_paths.links.overlap_nums[i] : 0) +
                gfaz_paths.links.overlap_ops[i];
    }
    for (auto& listener : this->edge_listeners) {
      listener(gfaz_paths.links.from_ids[i],
               i < gfaz_paths.links.from_orients.size() ? gfaz_paths.links.from_orients[i] == '-' : false,
               gfaz_paths.links.to_ids[i],
               i < gfaz_paths.links.to_orients.size() ? gfaz_paths.links.to_orients[i] == '-' : false,
               as_chars(overlap),
               GFAParser::tag_list_t());
    }
  }

  gfaz_paths.path_names =
      decompress_string_column(compressed.names_zstd, compressed.name_lengths_zstd);
  gfaz_paths.walk_sample_ids = decompress_string_column(
      compressed.walk_sample_ids_zstd, compressed.walk_sample_id_lengths_zstd);
  gfaz_paths.walk_hap_indices =
      Codec::zstd_decompress_uint32_vector(compressed.walk_hap_indices_zstd);
  gfaz_paths.walk_seq_ids = decompress_string_column(
      compressed.walk_seq_ids_zstd, compressed.walk_seq_id_lengths_zstd);
  gfaz_paths.walk_seq_starts = Codec::decompress_varint_int64(
      compressed.walk_seq_starts_zstd, compressed.walk_lengths.size());
  gfaz_paths.walk_seq_ends = Codec::decompress_varint_int64(
      compressed.walk_seq_ends_zstd, compressed.walk_lengths.size());

  auto decoded_rules = decode_rules(compressed);
  gfaz_paths.rules_first = std::move(decoded_rules.first);
  gfaz_paths.rules_second = std::move(decoded_rules.second);
  gfaz_paths.min_rule_id = compressed.min_rule_id();

  if (!compressed.paths_zstd.payload.empty()) {
    gfaz_paths.paths_flat = Codec::zstd_decompress_int32_vector(compressed.paths_zstd);
  }
  if (!compressed.walks_zstd.payload.empty()) {
    gfaz_paths.walks_flat = Codec::zstd_decompress_int32_vector(compressed.walks_zstd);
  }
  gfaz_paths.path_offsets = build_offsets(compressed.sequence_lengths);
  gfaz_paths.walk_offsets = build_offsets(compressed.walk_lengths);
  gfaz_paths.original_path_offsets = build_offsets(compressed.original_path_lengths);
  gfaz_paths.original_walk_offsets = build_offsets(compressed.original_walk_lengths);

  size_t encoded_path_count =
      gfaz_paths.path_offsets.empty() ? 0 : gfaz_paths.path_offsets.size() - 1;
  size_t path_count = min(gfaz_paths.path_names.size(), encoded_path_count);
  for (size_t i = 0; i < path_count; ++i) {
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.paths_flat, gfaz_paths.path_offsets,
        gfaz_paths.original_path_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, compressed.delta_round);
    auto visit_iteratee = make_gfaz_visit_iteratee(visits);
    auto overlap_iteratee = GFAParser::overlap_iteratee_t(
        [](const GFAParser::overlap_step_t&) {});
    for (auto& listener : this->path_listeners) {
      listener(gfaz_paths.path_names[i], visit_iteratee, overlap_iteratee, GFAParser::tag_list_t());
    }
  }

  size_t walk_count = gfaz_paths.walk_offsets.empty() ? 0 : gfaz_paths.walk_offsets.size() - 1;
  for (size_t i = 0; i < walk_count; ++i) {
    vector<NodeId> visits = decode_sequence_at_index(
        gfaz_paths.walks_flat, gfaz_paths.walk_offsets,
        gfaz_paths.original_walk_offsets, i, gfaz_paths.rules_first,
        gfaz_paths.rules_second, gfaz_paths.min_rule_id, compressed.delta_round);
    auto visit_iteratee = make_gfaz_visit_iteratee(visits);
    string sample_name = i < gfaz_paths.walk_sample_ids.size() ? gfaz_paths.walk_sample_ids[i] : "*";
    int64_t haplotype = i < gfaz_paths.walk_hap_indices.size() ? gfaz_paths.walk_hap_indices[i] : 0;
    string contig_name = i < gfaz_paths.walk_seq_ids.size() ? gfaz_paths.walk_seq_ids[i] : PathMetadata::NO_LOCUS_NAME;
    int64_t start = i < gfaz_paths.walk_seq_starts.size() ? gfaz_paths.walk_seq_starts[i] : PathMetadata::NO_END_POSITION;
    int64_t end = i < gfaz_paths.walk_seq_ends.size() ? gfaz_paths.walk_seq_ends[i] : PathMetadata::NO_END_POSITION;
    subrange_t subrange = PathMetadata::NO_SUBRANGE;
    if (start != PathMetadata::NO_END_POSITION) {
      subrange.first = start;
      if (end != PathMetadata::NO_END_POSITION) {
        subrange.second = end;
      }
    }
    for (auto& listener : this->walk_listeners) {
      listener(sample_name, haplotype, contig_name, subrange, visit_iteratee, GFAParser::tag_list_t());
    }
  }

  dispatch_rgfa_visits(gfaz_paths.node_count, gfaz_paths.segment_optional_fields,
                       gfaz_paths.node_lengths, this->max_rgfa_rank,
                       [&](nid_t id, int64_t offset, size_t length, const string& path_name, int64_t path_rank) {
                         for (auto& listener : this->rgfa_listeners) {
                           listener(id, offset, length, path_name, path_rank);
                         }
                       });
}

} // namespace algorithms
} // namespace vg
