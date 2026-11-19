/**
 * \file register_loader_saver_gfaz.cpp
 */

#include <vg/io/registry.hpp>
#include "register_loader_saver_gfaz.hpp"
#include "algorithms/gfa_to_handle.hpp"

#include "gfa.hpp"

#include "handle.hpp"
#include "save_handle_graph.hpp"
#include "../utility.hpp"

#include <GFAz/handlegraph_adapter/gfaz_handle_graph.hpp>

#include <cstdlib>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <unistd.h>

namespace vg {
namespace io {

using namespace std;
using namespace vg::io;

static bool sniff_gfaz(istream& stream) {
    if (!stream) {
        return false;
    }
    char buffer[sizeof(uint32_t)];
    size_t buffer_used = 0;
    while (stream.peek() != EOF && buffer_used < sizeof(buffer)) {
        buffer[buffer_used++] = (char) stream.get();
    }
    for (size_t i = 0; i < buffer_used; ++i) {
        stream.unget();
        if (!stream) {
            throw runtime_error("Ungetting failed after " + to_string(i) + " characters");
        }
    }
    if (!stream) {
        stream.clear();
    }
    return algorithms::GFAzParser::has_magic(buffer, buffer_used);
}

/// Spool a stream to a temporary file when the loader was not given a real
/// filename. Returns the name to read, and sets temp_name if one was made.
static string materialize_input(istream& input, const string& filename, string& temp_name) {
    if (!filename.empty() && filename != "-") {
        return filename;
    }
    temp_name = temp_file::create("gfaz-load");
    ofstream temp_stream(temp_name, ios::binary);
    temp_stream << input.rdbuf();
    temp_stream.close();
    return temp_name;
}

void register_loader_saver_gfaz() {
    // Read-only graph types are served straight from the container: paths stay
    // grammar-compressed and are decoded a block at a time. Set
    // VG_GFAZ_CONVERT=1 to get the old behavior of converting to a
    // GFAHandleGraph for every type.
    const char* convert_env = getenv("VG_GFAZ_CONVERT");
    bool always_convert = convert_env != nullptr && string(convert_env) != "0";

    if (!always_convert) {
        Registry::register_bare_loader_saver_with_header_check<gfaz::GfazHandleGraph, PathPositionHandleGraph, PathHandleGraph, HandleGraph>(
            "GFAZView", sniff_gfaz, [](istream& input, const string& filename) -> void* {
                string temp_name;
                try {
                    string load_name = materialize_input(input, filename, temp_name);
                    gfaz::GfazHandleGraphOptions options;
                    // Inherit vg's thread count instead of GFAz's own default.
                    options.num_threads = -1;
                    gfaz::GfazHandleGraph* graph = gfaz::GfazHandleGraph::load(load_name, options).release();
                    if (!temp_name.empty()) {
                        unlink(temp_name.c_str());
                    }
                    return (void*) graph;
                } catch (std::exception& e) {
                    cerr << "error[register_loader_saver_gfaz] GFAZ file " << filename
                         << " cannot be loaded: " << e.what() << endl;
                    exit(1);
                }
            }, [](const void* graph_void, ostream& output) {
                // The view is read-only and has no serialization of its own.
                cerr << "error[register_loader_saver_gfaz] a graph read from a GFAZ file "
                     << "cannot be saved directly; use vg convert" << endl;
                exit(1);
            });
    }

    auto convert_loader = [](istream& input, const string& filename) -> void* {
        GFAHandleGraph* gfa_graph = new GFAHandleGraph();
        string temp_name;
        try {
            string load_name = materialize_input(input, filename, temp_name);
            algorithms::GFAzParser parser;
            parser.external_id_map = &gfa_graph->gfa_id_space;
            algorithms::parser_to_path_handle_graph(parser, gfa_graph);
            parser.parse(load_name);
            gfa_graph->gfa_id_space.invert_translation();
        } catch (algorithms::GFAFormatError& e) {
            cerr << "error[register_loader_saver_gfaz] GFAZ ";
            if (!filename.empty() && filename != "-") {
                cerr << "file " << filename;
            } else {
                cerr << "stream";
            }
            cerr << " cannot be loaded: " << e.what() << endl;
            exit(1);
        }
        if (!temp_name.empty()) {
            unlink(temp_name.c_str());
        }
        return (void*) gfa_graph;
    };
    auto gfa_saver = [](const void* gfa_graph_void, ostream& output) {
        assert(gfa_graph_void != nullptr);
        graph_to_gfa((const GFAHandleGraph*) gfa_graph_void, output);
    };

    if (always_convert) {
        Registry::register_bare_loader_saver_with_header_check<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>(
            "GFAZ", sniff_gfaz, convert_loader, gfa_saver);
    } else {
        // Types that must be modified still need a real mutable graph.
        Registry::register_bare_loader_saver_with_header_check<GFAHandleGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph>(
            "GFAZ", sniff_gfaz, convert_loader, gfa_saver);
    }
}

}
}
