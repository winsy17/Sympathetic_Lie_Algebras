import unittest

class structure_constants_test(unittest.TestCase):
    def test_structure_constants(self):
        '''tests: 
        1) sl_2(C)

        The following tests come from the irreducible Lie algebras described in 
        the article 'Semi-invariants of low-dimensional Lie algebras' 
        https://doi.org/10.1080/00927872.2020.1861619
        We follow the same notation from the article.

        2) A_5_4
        3) A_5_5
        4) A_5_40 '''

        ##### 1 ##### 
        # sl_2(C)
        # Checking adjoints
        str_const_sl_2 = [[0,2,0],[0,0,-2],[1,0,0]]
        expected_adjoints_sl_2 = [[[0, 0, 0], [0, 2, 0], [0, 0, -2]], [[0, 0, 1], [-2, 0, 0], [0, 0, 0]], [[0, -1, 0], [0, 0, 0], [2, 0, 0]]]
        self.assertEqual(get_adjoints(str_const_sl_2),expected_adjoints_sl_2)
        # Checking is Lie algebra
        self.assertEqual(is_Lie_algebra(str_const_sl_2), bool(True))

        ##### 2 #####
        # A_5_4
        # Checking adjoints
        str_const_A_5_4 = [[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[1, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[1, 0, 0, 0, 0],[0, 0, 0, 0, 0]]
        expected_adjoints_A_5_4 = [[[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 1, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, -1, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, -1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]]

        self.assertEqual(get_adjoints(str_const_A_5_4),expected_adjoints_A_5_4)
        # Checking is Lie algebra
        self.assertEqual(is_Lie_algebra(str_const_A_5_4), bool(True))

        ##### 3 #####
        # A_5_5
        # Checking adjoints
        str_const_A_5_5 = [[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0],[1, 0, 0, 0, 0],[1, 0, 0, 0, 0],[0, 1, 0, 0, 0],[0, 0, 0, 0, 0]]
        expected_adjoints_A_5_5 = [[[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, -1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, -1, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]]

        self.assertEqual(get_adjoints(str_const_A_5_5),expected_adjoints_A_5_5)
        # Checking is Lie algebra
        self.assertEqual(is_Lie_algebra(str_const_A_5_5), bool(True))

        # ##### 4 #####
        # A_5_40
        # Checking adjoints
        str_const_A_5_40 = [[2, 0, 0, 0, 0],[0, -1, 0, 0, 0],[0, 0, 0, 0, 1],[0, 0, 0, 0, 0],[0, 0, 2, 0, 0],[0, 0, 0, 1, 0],[0, 0, 0, 0, -1],[0, 0, 0, 0, 0],[0, 0, 0, 1, 0],[0, 0, 0, 0, 0]]
        expected_adjoints_A_5_40 = [[[0, 2, 0, 0, 0], [0, 0, -1, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 1, 0]], [[-2, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 2, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, -1]], [[0, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, -2, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, -1, 0, 0, 0], [-1, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, -1, 0, 0], [0, 1, 0, 0, 0]]]

        self.assertEqual(get_adjoints(str_const_A_5_40),expected_adjoints_A_5_40)
        # Checking is Lie algebra
        self.assertEqual(is_Lie_algebra(str_const_A_5_40), bool(True))



if __name__=='__main__':
    from structure_constants import *
    unittest.main()