import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ecoli_cds_finder.cds_finder import get_coding_sequences


class TestCDSFinder(unittest.TestCase):
    def test_get_coding_sequences(self):
        test_seq = SeqRecord(Seq("AAAATGCCCCCTAGAAATGTTTTTTGAA"))
        expected_cds = ["ATGCCCCCTAG", "ATGTTTTTTGAA"]
        result = get_coding_sequences(test_seq)
        self.assertEqual(result, expected_cds)


if __name__ == "__main__":
    unittest.main()
