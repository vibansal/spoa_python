import unittest
import spoapy


def read_fastq(filename):
    res = []
    with open(filename) as inp:
        for line in inp:
            seq = next(inp).strip()
            _ = next(inp)
            qual = [ord(v) - 33 for v in next(inp).strip()]
            assert len(seq) == len(qual)
            res.append((seq, qual))
    return res


def allocate(engine, seqs):
    max_len = max(len(seq) for seq, _ in seqs)
    engine.preallocate(max_len, 4)


def update_graph(graph, engine, seqs, use_qualities):
    for seq, qual in seqs:
        if use_qualities:
            graph.add_sequence(engine, seq, qual)
        else:
            graph.add_sequence(engine, seq)


def calc_consensus(filename, use_quals, *args):
    seqs = read_fastq('test/data/sample.fastq')
    engine = spoapy.AlignmentEngine(*args)
    allocate(engine, seqs)

    graph = spoapy.Graph()
    update_graph(graph, engine, seqs, use_quals)
    consensus = graph.generate_consensus()
    return consensus


def test_msa(test, filename, use_quals, *args):
    seqs = read_fastq('test/data/sample.fastq')
    engine = spoapy.AlignmentEngine(*args)
    allocate(engine, seqs)

    graph = spoapy.Graph()
    update_graph(graph, engine, seqs, use_quals)
    msa = graph.generate_msa()
    test.assertEqual(len(msa), len(seqs))
    for msa_seq, (seq, _) in zip(msa, seqs):
        test.assertEqual(msa_seq.decode('utf-8').replace('-', ''), seq)


class TestSpoapy(unittest.TestCase):
    def test_graph_clear(self):
        seqs = read_fastq('test/data/sample.fastq')
        engine = spoapy.AlignmentEngine(spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -8, -8, -8)
        allocate(engine, seqs)

        graph = spoapy.Graph()
        update_graph(graph, engine, seqs, False)
        consensus1 = graph.generate_consensus()
        graph.clear()

        update_graph(graph, engine, seqs, False)
        consensus2 = graph.generate_consensus()
        self.assertEqual(consensus1, consensus2)

    def test_local_consensus(self):
        consensus = calc_consensus('test/data/sample.fastq', False,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -8, -8, -8)
        valid_result = "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGA"\
            "CCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGG"\
            "GAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGC"\
            "AGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGT"\
            "ACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGC"\
            "ACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTT"\
            "GAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATA"\
            "CGCTTTACACGCGCAACCAAGGATTTCGG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    def test_local_affine_consensus(self):
        consensus = calc_consensus('test/data/sample.fastq', False,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -6, -8, -6)
        valid_result = "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGA"\
            "CCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGG"\
            "GAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGC"\
            "AGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGT"\
            "ACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGC"\
            "ACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTT"\
            "GAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATA"\
            "CGCTTTACACGCGCAACCAAGGATTTCGG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    def test_local_convex_consensus(self):
        consensus = calc_consensus('test/data/sample.fastq', False,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -6, -10, -2)
        valid_result = "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGA"\
            "CCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGG"\
            "GAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGC"\
            "AGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGT"\
            "ACTCTGACACCGACGAATTTTACCCAGTTGCAGGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCG"\
            "CACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGT"\
            "TGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAAGCAGA"\
            "TACGCTG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    def test_local_consensus_qual(self):
        consensus = calc_consensus('test/data/sample.fastq', True,
            spoapy.AlignmentType.SW, 5, -4, -8, -8, -8, -8)
        valid_result = "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGA"\
            "CCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGG"\
            "GAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGC"\
            "AGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGT"\
            "ACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGC"\
            "ACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTT"\
            "GAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATA"\
            "CGCTTTACACGCGCAACCAAGGATTTCGG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    # Skipped two tests

    def test_global_consensus(self):
        consensus = calc_consensus('test/data/sample.fastq', False,
            spoapy.AlignmentType.NeedlemanWunsch, 5, -4, -8, -8, -8, -8)
        valid_result = "ATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGAC"\
            "CTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGG"\
            "AGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCA"\
            "GGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTA"\
            "CTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCA"\
            "CAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTG"\
            "AGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATAC"\
            "GC".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    def test_global_consensus_qual(self):
        consensus = calc_consensus('test/data/sample.fastq', True,
            spoapy.AlignmentType.NeedlemanWunsch, 5, -4, -8, -8, -8, -8)
        valid_result = "ATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGAC"\
            "CTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGG"\
            "AGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCA"\
            "GGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTA"\
            "CTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCA"\
            "CAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTG"\
            "AGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATAC"\
            "GC".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    # Skipped four tests

    def test_semiglobal_consensus(self):
        consensus = calc_consensus('test/data/sample.fastq', False,
            spoapy.AlignmentType.Overlap, 5, -4, -8, -8, -8, -8)
        valid_result = "ACATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCG"\
            "ACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAG"\
            "GGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGG"\
            "CAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCG"\
            "TACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCG"\
            "CACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGT"\
            "TGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGAT"\
            "ACGCGTTTTACACGCGCAACCAAGGATTTCGG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    def test_semiglobal_consensus_qual(self):
        consensus = calc_consensus('test/data/sample.fastq', True,
            spoapy.AlignmentType.Overlap, 5, -4, -8, -8, -8, -8)
        valid_result = "ACATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCG"\
            "ACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAG"\
            "GGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGG"\
            "CAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCG"\
            "TACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCG"\
            "CACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGT"\
            "TGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGAT"\
            "ACGCGTTTTACACGCGCAACCAAGGATTTCGG".encode('utf-8')
        self.assertEqual(consensus, valid_result)

    # Skipped four tests

    def test_msa(self):
        test_msa(self, 'test/data/sample.fastq', False,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -8, -8, -8)
        test_msa(self, 'test/data/sample.fastq', True,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -8, -8, -8)
        test_msa(self, 'test/data/sample.fastq', False,
            spoapy.AlignmentType.SmithWaterman, 5, -4, -8, -6, -8, -6)
        test_msa(self, 'test/data/sample.fastq', False,
            spoapy.AlignmentType.NeedlemanWunsch, 5, -4, -8, -8, -8, -8)
        test_msa(self, 'test/data/sample.fastq', False,
            spoapy.AlignmentType.Overlap, 5, -4, -8, -8, -8, -8)



if __name__ == "__main__":
    unittest.main()
