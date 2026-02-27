#ifndef VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED
#define VG_ALGORITHMS_GFAZ_TO_HANDLE_HPP_INCLUDED

/**
 * \file gfaz_to_handle.hpp
 *
 * Defines algorithms for copying data from GFAZ files into handle graphs.
 */

#include "gfa_to_handle.hpp"

namespace vg {
namespace algorithms {
using namespace std;

/// Read a GFAZ file into a HandleGraph.
void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          GFAIDMapInfo* translation = nullptr,
                          int num_threads = 0);

/// Overload which serializes its translation to a file internally.
void gfaz_to_handle_graph(const string& filename,
                          MutableHandleGraph* graph,
                          const string& translation_filename,
                          int num_threads = 0);

/// Read a GFAZ file into a PathHandleGraph.
void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               GFAIDMapInfo* translation = nullptr,
                               int64_t max_rgfa_rank = numeric_limits<int64_t>::max(),
                               unordered_set<PathSense>* ignore_sense = nullptr,
                               int num_threads = 0);

/// Overload which serializes its translation to a file internally.
void gfaz_to_path_handle_graph(const string& filename,
                               MutablePathMutableHandleGraph* graph,
                               int64_t max_rgfa_rank,
                               const string& translation_filename,
                               unordered_set<PathSense>* ignore_sense = nullptr,
                               int num_threads = 0);

}
}

#endif
