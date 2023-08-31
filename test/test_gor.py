import unittest

from source.GOR import *
from source.Parser import *

directory = "ressources/dataset/dssp"


class TestGOR(unittest.TestCase):
    def setUp(self):
        self.parser = DSSPParser(directory)
        self.gor = GOR(self.parser)

    def test_MCC(self):
        self.assertEqual(self.gor.MCC("EEEEEEHHCEEE", "CCHEEEHECEEC", "H"), 40.0,  msg = "Check how you count the TP/TN/FP/FN, it depends on the structure !")
        self.assertEqual(self.gor.MCC("EEEEEEHHCEEE", "CCHEEEHECEEC", "E"), 19.25)
        self.assertEqual(self.gor.MCC("EEEEEEHHCEEE", "CCHEEEHECEEC", "C"), 42.64)

    def test_overall_performance(self):
        Q3,MCCH,MCCE,MCCC = self.gor.validate(self.parser.get_testset())
        self.assertTrue(60 < Q3[0] < 63, msg = "Check the counters, did you take into account the position of the"
                                            "amino acid in the neighbourhood ? This is the most common mistake")
        self.assertTrue(35 < MCCH[0] < 41, msg = "Check how you count the TP/TN/FP/FN, it depends on the structure !")
        self.assertTrue(35 < MCCE[0] < 41, msg = "Don't forget your mean should exclude the division errors")
        self.assertTrue(35 < MCCC[0] < 41)



if __name__ == '__main__':
    unittest.main()
