from libc.stdint cimport int8_t, uint8_t, int32_t, uint32_t, int64_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr
from libcpp.pair cimport pair
from libcpp.string cimport string


ctypedef vector[pair[int32_t, int32_t]] Alignment

cdef extern from 'graph.hpp' namespace 'spoa':
    unique_ptr[Graph] createGraph()

    cdef cppclass Node:
        uint32_t id()
        uint32_t code()
        vector[shared_ptr[Edge]]& in_edges()
        vector[shared_ptr[Edge]]& out_edges()
        vector[uint32_t]& aligned_nodes_ids()
        bool successor(uint32_t& dst, uint32_t label)
        uint32_t coverage()

    cdef cppclass Edge:
        uint32_t begin_node_id()
        uint32_t end_node_id()
        int64_t weight()

    cdef cppclass Graph:
        vector[unique_ptr[Node]]& nodes()
        vector[uint32_t]& rank_to_node_id()
        uint32_t num_codes()
        uint8_t coder(uint8_t c)
        uint8_t decoder(uint8_t code)

        void add_alignment(Alignment& alignment, string& sequence, uint32_t weight)
        void add_alignment(Alignment& alignment, string& sequence, vector[uint32_t]& weights)

        void generate_multiple_sequence_alignment(vector[string]& dst, bool include_consensus)
        string generate_consensus()
        string generate_consensus(vector[uint32_t]& dst, bool verbose)

        unique_ptr[Graph] subgraph(uint32_t begin_node_id, uint32_t end_node_id,
            vector[int32_t]& subgraph_to_graph_mapping)
        void update_alignment(Alignment& alignment, vector[int32_t]& subgraph_to_graph_mapping)
        void print_dot(string& path)
        vector[shared_ptr[Edge]] consensus_edges()
        void clear()


cdef extern from "<utility>" namespace "std" nogil:
    cdef unique_ptr[Graph] move(unique_ptr[Graph])
    cdef unique_ptr[Alignment] move(unique_ptr[Alignment])
    cdef Alignment move(Alignment)


cdef extern from 'alignment_engine.hpp' namespace 'spoa':
    cdef enum _AlignmentType 'spoa::AlignmentType':
        _kSW 'spoa::AlignmentType::kSW'
        _kNW 'spoa::AlignmentType::kNW'
        _kOV 'spoa::AlignmentType::kOV'

    cpdef enum _AlignmentSubtype 'spoa::AlignmentSubtype':
        _kLinear 'spoa::AlignmentSubtype::kLinear',
        _kAffine 'spoa::AlignmentSubtype::kAffine',
        _kConvex 'spoa::AlignmentSubtype::kConvex'

    cdef cppclass AlignmentEngine:
        void prealloc(uint32_t max_sequence_size, uint32_t alphabet_size)
        Alignment align(string& sequence, unique_ptr[Graph]& graph)

    unique_ptr[AlignmentEngine] createAlignmentEngine(_AlignmentType type,
            int8_t m, int8_t n, int8_t g, int8_t e, int8_t q, int8_t c)

cpdef enum AlignmentType:
    kSW = <uint32_t> _kSW
    kNW = <uint32_t> _kNW
    kOV = <uint32_t> _kOV

cpdef enum AlignmentSubtype:
    kLinear = <uint32_t> _kLinear
    kAffine = <uint32_t> _kAffine
    kConvex = <uint32_t> _kConvex
