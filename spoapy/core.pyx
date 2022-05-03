# distutils: language = c++
# distutils: sources = src/graph.cpp src/alignment_engine.cpp

from spoapy cimport cspoa
from libc.stdint cimport int8_t, uint8_t, int32_t, uint32_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr, make_unique
from libcpp.pair cimport pair
from libcpp.string cimport string
from cython.operator cimport dereference
from enum import Enum

cpdef enum AlignmentType:
    SmithWaterman = cspoa.AlignmentType.kSW
    NeedlemanWunsch = cspoa.AlignmentType.kNW
    Overlap = cspoa.AlignmentType.kOV
    SW = cspoa.AlignmentType.kSW
    NW = cspoa.AlignmentType.kNW
    OV = cspoa.AlignmentType.kOV


cpdef enum AlignmentSubtype:
    Linear = cspoa.AlignmentSubtype.kLinear
    Affine = cspoa.AlignmentSubtype.kAffine
    Convex = cspoa.AlignmentSubtype.kConvex


def seq_to_cstr(sequence):
    if isinstance(sequence, bytes):
        return <string>sequence
    else:
        return <string>sequence.encode('utf-8')


cdef class AlignmentEngine:
    cdef unique_ptr[cspoa.AlignmentEngine] _c_aln_engine

    def __cinit__(self, AlignmentType type,
            int8_t m=5, int8_t n=-4, int8_t g=-8, int8_t e=-6, int8_t q=-10, int8_t c=-4):
        """
        Create an alignment engine:

        Parameters
        ----------
        type: AlignmentType
        m: int
            default: 5
            score for matching bases
        n: int
            default: -4
            score for mismatching bases
        g: int
            default: -8
            gap opening penalty (must be non-positive)
        e: int
            default: -6
            gap extension penalty (must be non-positive)
        q: int
            default: -10
            gap opening penalty of the second affine function
            (must be non-positive)
        c: int
            default: -4
            gap extension penalty of the second affine function
            (must be non-positive)
        """
        self._c_aln_engine = cspoa.createAlignmentEngine(<cspoa._AlignmentType>type, m, n, g, e, q, c)

    def preallocate(self, uint32_t max_sequence_size, uint32_t alphabet_size):
        self._c_aln_engine.get().prealloc(max_sequence_size, alphabet_size)


cdef class Graph:
    cdef unique_ptr[cspoa.Graph] _c_graph

    @staticmethod
    cdef _init(unique_ptr[cspoa.Graph] graph):
        cdef Graph self = Graph.__new__(Graph)
        self._c_graph = cspoa.move(graph)
        return self

    def __len__(self):
        """
        Returns the number of nodes.
        """
        return self._c_graph.get().nodes().size()

    def __cinit__(self):
        self._c_graph = cspoa.createGraph()

    def nodes(self):
        return NodeVector._init(&self._c_graph.get().nodes())

    def rank_to_node_id(self):
        return IntVector._init(&self._c_graph.get().rank_to_node_id())

    def alphabet_size(self):
        """
        Returns the number of unique letters in the alignment graph.
        Named `num_codes` in the spoa library.
        """
        return self._c_graph.get().num_codes()

    def coder(self, uint8_t c):
        return self._c_graph.get().coder(c)

    def decoder(self, uint8_t code):
        return self._c_graph.get().decoder(code)

    def update_edge(self, start_node_id, end_node_id, weight):
            return self._c_graph.get().update_edge(<uint32_t>start_node_id, <uint32_t>end_node_id, <int32_t>weight)

    def get_path_weights(self, n, nodes):
         cdef vector[uint32_t] nodes_v = nodes
         return self._c_graph.get().get_path_weights(<uint32_t>n, nodes_v)

    def align(self, AlignmentEngine engine, sequence):
        cdef unique_ptr[cspoa.Alignment] aln_ptr = unique_ptr[cspoa.Alignment](
            new cspoa.Alignment(engine._c_aln_engine.get().align(seq_to_cstr(sequence), self._c_graph)))
        return Alignment._init(cspoa.move(aln_ptr))

    def add_alignment(self, Alignment alignment, sequence, weight=1):
        """
        Adds an alignment to the graph.

        Parameters
        ----------
        alignment: Alignment
        sequence: bytes or str
        weight: int or list[int]
            default: 1
            Weight of the sequence (single number) or weight of each individual letter (same length as the sequence).
        """
        cdef string cstr = seq_to_cstr(sequence)
        try:
            assert len(weight) == len(sequence)
        except TypeError:
            return self._c_graph.get().add_alignment(dereference(alignment._c_vec.get()), cstr, <uint32_t>weight)

        cdef vector[uint32_t] weights = weight
        return self._c_graph.get().add_alignment(dereference(alignment._c_vec.get()), cstr, weights)

    def add_sequence(self, AlignmentEngine engine, sequence, weight=1):
        """
        Aligns the sequence and adds it to the graph. Same as calling `align` and `add_alignment`.

        Parameters
        ----------
        engine: AlignmentEngine
        sequence: bytes or str
        weight: int or list[int]
            default: 1
            Weight of the sequence (single number) or weight of each individual letter (same length as the sequence).
        """
        cdef string cstr = seq_to_cstr(sequence)
        cdef cspoa.Alignment alignment = engine._c_aln_engine.get().align(cstr, self._c_graph)

        try:
            assert len(weight) == len(sequence)
        except TypeError:
            return self._c_graph.get().add_alignment(alignment, cstr, <uint32_t>weight)

        cdef vector[uint32_t] weights = weight
        return self._c_graph.get().add_alignment(alignment, cstr, weights)

    def generate_msa(self, bool include_consensus=False):
        """
        Generates multiple sequence alignment. Returns list[bytes].
        """

        cdef vector[string] dst
        self._c_graph.get().generate_multiple_sequence_alignment(dst, include_consensus)
        n = dst.size()
        return [dst.at(i) for i in range(n)]

    def generate_consensus(self, bool coverage=False, bool matrix=False):
        """
        Generates consensus.

        Returns:
        * single bytes string with consensus;
        * if `coverage` is true, returns tuple (bytes, list[int]) where the second element is the base coverage;
        * if `matrix` is true, returns tuple (bytes, list[list[int]]) where the second element is
            the complete summary matrix (alphabet_size + 1) x consensus_size.

        Note, that `coverage` and `matrix` cannot be both true.
        """

        assert not (coverage and matrix)
        if not coverage and not matrix:
            return self._c_graph.get().generate_consensus()

        cdef vector[uint32_t] dst
        consensus = self._c_graph.get().generate_consensus(dst, matrix)

        n = len(consensus)
        if coverage:
            assert dst.size() == n
            return consensus, dst

        m = self.alphabet_size() + 1
        assert dst.size() == m * n
        res = [[0] * n for _ in range(m)]
        for i in range(m):
            for j in range(n):
                res[i][j] = dst.at(i * n + j)
        return consensus, res

    # def subgraph(self, uint32_t start_node_id, uint32_t end_node_id):
    #     """
    #     Creates a subgraph and returns a tuple (Graph, list[int]),
    #     where the second element represents subgraph to graph mapping.
    #     """
    #     cdef vector[int32_t] subgraph_to_graph
    #     cdef Graph subgraph = Graph._init(self._c_graph.get().subgraph(start_node_id, end_node_id, subgraph_to_graph))
    #     return subgraph, subgraph_to_graph

    def print_dot(self, path):
        self._c_graph.get().print_dot(<string>path.encode('utf-8'))

    def consensus_edges(self):
        """
        Returns consensus and the a of edges.
        """
        consensus = self._c_graph.get().generate_consensus()
        cdef vector[shared_ptr[cspoa.Edge]] edges = self._c_graph.get().consensus_edges()
        return consensus, [Edge._init(edges.at(i)) for i in range(edges.size())]

    def clear(self):
        self._c_graph.get().clear()


def __pri_ctor(self):
    raise TypeError('Class "%s" cannot be constructed' % self.__class__.__name__)


cdef class Node:
    cdef unique_ptr[cspoa.Node]* _c_node

    @staticmethod
    cdef _init(unique_ptr[cspoa.Node]* node):
        cdef Node self = Node.__new__(Node)
        self._c_node = node
        return self

    def __init__(self):
        __pri_ctor(self)

    def id(self):
        return self._c_node.get().id()

    def code(self):
        return self._c_node.get().code()

    def in_edges(self):
        return EdgeVector._init(&self._c_node.get().in_edges())

    def out_edges(self):
        return EdgeVector._init(&self._c_node.get().out_edges())

    def aligned_nodes_ids(self):
        return IntVector._init(&self._c_node.get().aligned_nodes_ids())

    def successor(self, uint32_t label):
        cdef uint32_t dst = 0
        return dst if self._c_node.get().successor(dst, label) else None

    def coverage(self):
        return self._c_node.get().coverage()


cdef class Edge:
    cdef shared_ptr[cspoa.Edge] _c_edge

    @staticmethod
    cdef _init(shared_ptr[cspoa.Edge] edge):
        cdef Edge self = Edge.__new__(Edge)
        self._c_edge = edge
        return self

    def __init__(self):
        __pri_ctor(self)

    def start_node_id(self):
        return self._c_edge.get().begin_node_id()

    def end_node_id(self):
        return self._c_edge.get().end_node_id()

    def weight(self):
        return self._c_edge.get().weight()

    def labels(self):
        return self._c_edge.get().labels()


# All types of vectors below.


cdef class Alignment:
    cdef unique_ptr[cspoa.Alignment] _c_vec

    @staticmethod
    cdef _init(unique_ptr[cspoa.Alignment] vec):
        cdef Alignment self = Alignment.__new__(Alignment)
        self._c_vec = cspoa.move(vec)
        return self

    def __init__(self):
        __pri_ctor(self)

    def __len__(self):
        return self._c_vec.get().size()

    def __getitem__(self, size_t index):
        return self._c_vec.get().at(index)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    @staticmethod
    def from_pairs(pairs):
        cdef Alignment self = Alignment.__new__(Alignment)
        self._c_vec = unique_ptr[cspoa.Alignment](new cspoa.Alignment())
        for x, y in pairs:
            z = pair[int32_t, int32_t](x, y)
            self._c_vec.get().push_back(z)
        return self

    def __str__(self):
        res = 'Alignment['
        res += ', '.join('(%d,%d)' % (x, y) for x, y in self)
        res += ']'
        return res


cdef class IntVector:
    cdef const vector[uint32_t]* _c_vec

    @staticmethod
    cdef _init(const vector[uint32_t]* vec):
        cdef IntVector self = IntVector.__new__(IntVector)
        self._c_vec = vec
        return self

    def __init__(self):
        __pri_ctor(self)

    def __len__(self):
        return self._c_vec.size()

    def __getitem__(self, size_t index):
        return self._c_vec.at(index)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


cdef class NodeVector:
    cdef const vector[unique_ptr[cspoa.Node]]* _c_vec

    @staticmethod
    cdef _init(const vector[unique_ptr[cspoa.Node]]* vec):
        cdef NodeVector self = NodeVector.__new__(NodeVector)
        self._c_vec = vec
        return self

    def __init__(self):
        __pri_ctor(self)

    def __len__(self):
        return self._c_vec.size()

    def __getitem__(self, size_t index):
        return Node._init(&self._c_vec.at(index))

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


cdef class EdgeVector:
    cdef const vector[shared_ptr[cspoa.Edge]]* _c_vec

    @staticmethod
    cdef _init(const vector[shared_ptr[cspoa.Edge]]* vec):
        cdef EdgeVector self = EdgeVector.__new__(EdgeVector)
        self._c_vec = vec
        return self

    def __init__(self):
        __pri_ctor(self)

    def __len__(self):
        return self._c_vec.size()

    def __getitem__(self, size_t index):
        return Edge._init(self._c_vec.at(index))

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


__all__ = ['AlignmentEngine', 'Graph', 'AlignmentType', 'AlignmentSubtype', 'Alignment']
