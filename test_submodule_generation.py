import unittest

class submodule_generation_tests(unittest.TestCase):
    def test_gen_repr_tensor(self):
        '''tests: 
        1) generating the simple module V_4 in the tensor product of V_2 tensor V_2
        2) generating the simple module V_2 in the tensor product of V_2 tensor V_2
        2) generating the simple module V_0 in the tensor product of V_2 tensor V_2'''
        n = 2
        m = 2
        k = 4
        space = 'tensor'
        
        # Check generation of elements of weight 4 
        expected_weight_4 = [(0, 0)]
        self.assertEqual(get_fixed_weight(n,m,k,space),expected_weight_4)

        # Check generation of elements of weight 6
        expected_weight_6 = []
        self.assertEqual(get_fixed_weight(n,m,k+2,space),expected_weight_6)

        # Check generation of elements of weights 2,0,-2
        expected_weight_2 = [(0, 1), (1, 0)]
        expected_weight_0 = [(0, 2), (1, 1),(2, 0)]
        expected_weight__2 = [(1, 2), (2, 1)]
        self.assertEqual(get_fixed_weight(n,m,2,space),expected_weight_2)
        self.assertEqual(get_fixed_weight(n,m,0,space),expected_weight_0)
        self.assertEqual(get_fixed_weight(n,m,-2,space),expected_weight__2)

        # Check rho_E 
        # in elements of weight 4
        self.assertEqual(rho_E(n,m,space,(0,0)),[])
        # in elements of weight 2
        self.assertEqual(rho_E(n,m,space,(1,0)),[[2, 0, 0]])
        self.assertEqual(rho_E(n,m,space,(0,1)),[[2, 0, 0]])
        # in elements of weight 0
        self.assertEqual(rho_E(n,m,space,(2,0)),[[2, 1, 0]])
        self.assertEqual(rho_E(n,m,space,(1,1)),[[2, 0, 1], [2, 1, 0]])
        self.assertEqual(rho_E(n,m,space,(0,2)),[[2, 0, 1]])
        # in elements of weight -2
        self.assertEqual(rho_E(n,m,space,(1,2)),[[2, 0, 2], [2, 1, 1]])
        self.assertEqual(rho_E(n,m,space,(2,1)),[[2, 1, 1],[2, 2, 0]])

        # Check kernel of rho_E
        expected_ker_w_4 = [1]
        self.assertListEqual(kernel_rho_E(expected_weight_4,expected_weight_6,n,m,space),expected_ker_w_4)
        
        expected_ker_w_2 = [1, -1]
        self.assertListEqual(kernel_rho_E(expected_weight_2,expected_weight_4,n,m,space),expected_ker_w_2)

        # Check rho_F
        # weight 4
        expected_rho_F_00 = [[1, 1, 0], [1, 0, 1]]
        self.assertListEqual(rho_F(n,m,space,(0,0)),expected_rho_F_00)
        # weight 2
        expected_rho_F_10 = [[1, 2, 0], [1, 1, 1]]
        self.assertListEqual(rho_F(n,m,space,(1,0)),expected_rho_F_10)
        expected_rho_F_01 = [[1, 1, 1], [1, 0, 2]]
        self.assertListEqual(rho_F(n,m,space,(0,1)),expected_rho_F_01)
        # weight 0
        expected_rho_F_02 = [[1, 1, 2]]
        self.assertListEqual(rho_F(n,m,space,(0,2)),expected_rho_F_02)
        expected_rho_F_11 = [[1, 2, 1],[1, 1, 2]]
        self.assertListEqual(rho_F(n,m,space,(1,1)),expected_rho_F_11)
        expected_rho_F_20 = [[1, 2, 1]]
        self.assertListEqual(rho_F(n,m,space,(2,0)),expected_rho_F_20)

        # Check linear_rho_F
        u_4 = [[1,0,0]]
        expected_rho_w_4 = [[1, 1, 0], [1, 0, 1]]
        self.assertListEqual(linear_rho_F(u_4,n,m,space),expected_rho_w_4)
        u_2 = [[-1,0,1],[1,1,0]]
        expected_rho_w_2 = [[-1, 0, 2], [1, 2, 0]]
        self.assertListEqual(linear_rho_F(u_2,n,m,space),expected_rho_w_2)
        u_0 = [[1,0,2],[-1,1,1],[1,2,0]]
        expected_rho_w_0 = []
        self.assertListEqual(linear_rho_F(u_0,n,m,space),expected_rho_w_0)

        # Check rho_F_iterates
        expected_iterates_w_4 = [[[1, 0, 0]], [[1, 1, 0], [1, 0, 1]], [[1, 2, 0], [2, 1, 1], [1, 0, 2]], [[3, 2, 1], [3, 1, 2]], [[6, 2, 2]]]
        self.assertListEqual(rho_F_iterates(u_4,n,m,space,4),expected_iterates_w_4)

        expected_iterates_w_2 = [[[-1, 0, 1], [1, 1, 0]], [[-1, 0, 2], [1, 2, 0]], [[-1, 1, 2], [1, 2, 1]]]
        self.assertListEqual(rho_F_iterates(u_2,n,m,space,2),expected_iterates_w_2)

        expected_iterates_w_0 = [[[1, 0, 2], [-1, 1, 1], [1, 2, 0]]]
        self.assertListEqual(rho_F_iterates(u_0,n,m,space,0),expected_iterates_w_0)

        # Check gen_repr
        expected_repr_4 = [[[1, 0, 0]], [[1, 1, 0], [1, 0, 1]], [[1, 2, 0], [2, 1, 1], [1, 0, 2]], [[3, 2, 1], [3, 1, 2]], [[6, 2, 2]]]
        self.assertListEqual(gen_repr(n,m,k,space),expected_repr_4)

        expected_repr_2 = [[[1, 0, 1], [-1, 1, 0]], [[1, 0, 2], [-1, 2, 0]], [[1, 1, 2], [-1, 2, 1]]]
        self.assertListEqual(gen_repr(n,m,2,space),expected_repr_2)

        expected_repr_0 = [[[1.0, 0, 2], [-1.0, 1, 1], [1.0, 2, 0]]]
        self.assertListEqual(gen_repr(n,m,0,space),expected_repr_0)
    
    def test_gen_repr_symmetric(self):
        '''tests: 
        1) generating the simple modules V_12, V_8, V_4, V_0 in the symmetric power S^2(V_6)'''
        n = 6
        m = 6
        space = 'symmetric power'
        
        # Check generation of elements of weight 12,8,4,0 
        expected_weight_14 = []
        expected_weight_12 = [(0, 0)]
        expected_weight_8 = [(0, 2), (1, 1)]
        expected_weight_4 = [(0, 4), (1, 3), (2, 2)]
        expected_weight_0 = [(0, 6), (1, 5), (2, 4), (3, 3)]
        
        self.assertEqual(get_fixed_weight(n,m,12,space),expected_weight_12)

        self.assertEqual(get_fixed_weight(n,m,8,space),expected_weight_8)
        
        self.assertEqual(get_fixed_weight(n,m,4,space),expected_weight_4)

        self.assertEqual(get_fixed_weight(n,m,0,space),expected_weight_0)

        expected_weight_10 = [(0, 1)]
        self.assertEqual(get_fixed_weight(n,m,10,space),expected_weight_10)
        expected_weight_6 = [(0, 3), (1, 2)]
        self.assertEqual(get_fixed_weight(n,m,6,space),expected_weight_6)
        expected_weight_2 = [(0, 5), (1, 4), (2, 3)]
        self.assertEqual(get_fixed_weight(n,m,2,space),expected_weight_2)


        # Check rho_E 
        # in elements of weight 12
        self.assertEqual(rho_E(n,m,space,(0, 0)), [] )
        
        # in elements of weight 8
        self.assertEqual(rho_E(n,m,space,(0, 2)), [[10, 0, 1]] )
        self.assertEqual(rho_E(n,m,space,(1, 1)), [[6, 0, 1]] )
        
        # in elements of weight 4
        self.assertEqual(rho_E(n,m,space,(0, 4)), [[12, 0, 3]] )
        self.assertEqual(rho_E(n,m,space,(1, 3)), [[6, 0, 3], [12, 1, 2]] )
        self.assertEqual(rho_E(n,m,space,(2, 2)), [[10, 1, 2]] )
        
        # in elements of weight 0
        self.assertEqual(rho_E(n,m,space,(0, 6)), [[6, 0, 5]] )
        self.assertEqual(rho_E(n,m,space,(1, 5)), [[6, 0, 5], [10, 1, 4]] )
        self.assertEqual(rho_E(n,m,space,(2, 4)), [[10, 1, 4], [12, 2, 3]] )
        self.assertEqual(rho_E(n,m,space,(3, 3)), [[12, 2, 3]] )

        # Check kernel of rho_E
        expected_ker_w_12 =  [1]
        expected_ker_w_8 =  [3, -5]
        expected_ker_w_4 =  [5, -10, 12]
        expected_ker_w_0 =  [1, -1, 1, -1]

        self.assertListEqual(kernel_rho_E(expected_weight_12,expected_weight_14,n,m,space),expected_ker_w_12)
        self.assertListEqual(kernel_rho_E(expected_weight_8,expected_weight_10,n,m,space),expected_ker_w_8)
        self.assertListEqual(kernel_rho_E(expected_weight_4,expected_weight_6,n,m,space),expected_ker_w_4)
        self.assertListEqual(kernel_rho_E(expected_weight_0,expected_weight_2,n,m,space),expected_ker_w_0)


        # Check rho_F
        # weight  12
        expected_rho_F_00 =  [[1, 0, 1]]
        self.assertListEqual(rho_F(n,m,space,(0, 0)),expected_rho_F_00)
        # weight  8
        expected_rho_F_02 =  [[1, 1, 2], [1, 0, 3]]
        self.assertListEqual(rho_F(n,m,space,(0, 2)),expected_rho_F_02)
        expected_rho_F_11 =  [[1, 1, 2]]
        self.assertListEqual(rho_F(n,m,space,(1, 1)),expected_rho_F_11)
        # weight  4
        expected_rho_F_04 =  [[1, 1, 4], [1, 0, 5]]
        self.assertListEqual(rho_F(n,m,space,(0, 4)),expected_rho_F_04)
        expected_rho_F_13 =  [[1, 2, 3], [1, 1, 4]]
        self.assertListEqual(rho_F(n,m,space,(1, 3)),expected_rho_F_13)
        expected_rho_F_22 =  [[1, 2, 3]]
        self.assertListEqual(rho_F(n,m,space,(2, 2)),expected_rho_F_22)
        # weight  0
        expected_rho_F_06 =  [[1, 1, 6]]
        self.assertListEqual(rho_F(n,m,space,(0, 6)),expected_rho_F_06)
        expected_rho_F_15 =  [[1, 2, 5], [1, 1, 6]]
        self.assertListEqual(rho_F(n,m,space,(1, 5)),expected_rho_F_15)
        expected_rho_F_24 =  [[1, 3, 4], [1, 2, 5]]
        self.assertListEqual(rho_F(n,m,space,(2, 4)),expected_rho_F_24)
        expected_rho_F_33 =  [[1, 3, 4]]
        self.assertListEqual(rho_F(n,m,space,(3, 3)),expected_rho_F_33)

        # Check linear_rho_F
        # weight  12
        u_12 =  [[1, 0, 0]]
        expected_rho_w_12 =  [[1, 0, 1]]
        self.assertListEqual(linear_rho_F(u_12,n,m,space),expected_rho_w_12)
        # weight  8
        u_8 =  [[3, 0, 2], [-5, 1, 1]]
        expected_rho_w_8 =  [[-2, 1, 2], [3, 0, 3]]
        self.assertListEqual(linear_rho_F(u_8,n,m,space),expected_rho_w_8)
        # weight  4
        u_4 =  [[5, 0, 4], [-10, 1, 3], [12, 2, 2]]
        expected_rho_w_4 =  [[-5, 1, 4], [5, 0, 5], [2, 2, 3]]
        self.assertListEqual(linear_rho_F(u_4,n,m,space),expected_rho_w_4)
        # weight  0
        u_0 =  [[1, 0, 6], [-1, 1, 5], [1, 2, 4], [-1, 3, 3]]
        expected_rho_w_0 =  []
        self.assertListEqual(linear_rho_F(u_0,n,m,space),expected_rho_w_0)

        # Check rho_F_iterates
        expected_iterates_w_12 =  [[[1, 0, 0]], [[1, 0, 1]], [[2, 1, 1], [1, 0, 2]], [[3, 1, 2], [1, 0, 3]], [[6, 2, 2], [4, 1, 3], [1, 0, 4]], [[10, 2, 3], [5, 1, 4], [1, 0, 5]], [[20, 3, 3], [15, 2, 4], [6, 1, 5], [1, 0, 6]], [[35, 3, 4], [21, 2, 5], [7, 1, 6]], [[70, 4, 4], [56, 3, 5], [28, 2, 6]], [[126, 4, 5], [84, 3, 6]], [[252, 5, 5], [210, 4, 6]], [[462, 5, 6]], [[924, 6, 6]]]
        self.assertListEqual(rho_F_iterates(u_12,n,m,space,12),expected_iterates_w_12)
        expected_iterates_w_8 =  [[[3, 0, 2], [-5, 1, 1]], [[-2, 1, 2], [3, 0, 3]], [[-4, 2, 2], [1, 1, 3], [3, 0, 4]], [[-3, 2, 3], [4, 1, 4], [3, 0, 5]], [[-6, 3, 3], [1, 2, 4], [7, 1, 5], [3, 0, 6]], [[-5, 3, 4], [8, 2, 5], [10, 1, 6]], [[-10, 4, 4], [3, 3, 5], [18, 2, 6]], [[-7, 4, 5], [21, 3, 6]], [[-14, 5, 5], [14, 4, 6]]]
        self.assertListEqual(rho_F_iterates(u_8,n,m,space,8),expected_iterates_w_8)
        expected_iterates_w_4 =  [[[5, 0, 4], [-10, 1, 3], [12, 2, 2]], [[-5, 1, 4], [5, 0, 5], [2, 2, 3]], [[-3, 2, 4], [5, 0, 6], [4, 3, 3]], [[1, 3, 4], [-3, 2, 5], [5, 1, 6]], [[2, 4, 4], [-2, 3, 5], [2, 2, 6]]]
        self.assertListEqual(rho_F_iterates(u_4,n,m,space,4),expected_iterates_w_4)
        expected_iterates_w_0 =  [[[1, 0, 6], [-1, 1, 5], [1, 2, 4], [-1, 3, 3]]]
        self.assertListEqual(rho_F_iterates(u_0,n,m,space,0),expected_iterates_w_0)

        # Check gen_repr
        expected_repr_12 = [[[1, 0, 0]], [[1, 0, 1]], [[2, 1, 1], [1, 0, 2]], [[3, 1, 2], [1, 0, 3]], [[6, 2, 2], [4, 1, 3], [1, 0, 4]], [[10, 2, 3], [5, 1, 4], [1, 0, 5]], [[20, 3, 3], [15, 2, 4], [6, 1, 5], [1, 0, 6]], [[35, 3, 4], [21, 2, 5], [7, 1, 6]], [[70, 4, 4], [56, 3, 5], [28, 2, 6]], [[126, 4, 5], [84, 3, 6]], [[252, 5, 5], [210, 4, 6]], [[462, 5, 6]], [[924, 6, 6]]]
        self.assertListEqual(gen_repr(n,m,12,space), expected_repr_12)
        expected_repr_8 = [[[3, 0, 2], [-5, 1, 1]], [[-2, 1, 2], [3, 0, 3]], [[-4, 2, 2], [1, 1, 3], [3, 0, 4]], [[-3, 2, 3], [4, 1, 4], [3, 0, 5]], [[-6, 3, 3], [1, 2, 4], [7, 1, 5], [3, 0, 6]], [[-5, 3, 4], [8, 2, 5], [10, 1, 6]], [[-10, 4, 4], [3, 3, 5], [18, 2, 6]], [[-7, 4, 5], [21, 3, 6]], [[-14, 5, 5], [14, 4, 6]]]
        self.assertListEqual(gen_repr(n,m,8,space), expected_repr_8)
        expected_repr_4 = [[[5, 0, 4], [-10, 1, 3], [12, 2, 2]], [[-5, 1, 4], [5, 0, 5], [2, 2, 3]], [[-3, 2, 4], [5, 0, 6], [4, 3, 3]], [[1, 3, 4], [-3, 2, 5], [5, 1, 6]], [[2, 4, 4], [-2, 3, 5], [2, 2, 6]]]
        self.assertListEqual(gen_repr(n,m,4,space), expected_repr_4)
        expected_repr_0 = [[[1, 0, 6], [-1, 1, 5], [1, 2, 4], [-1, 3, 3]]]
        self.assertListEqual(gen_repr(n,m,0,space), expected_repr_0)

    def test_gen_repr_exterior(self):
        '''tests: 
        generating the simple modules V_10, V_6, V_2 in the exterior power \Lambda^2(V_6)'''
        n = 6
        m = 6
        space = 'exterior power'
        
        # Check generation of elements of weight 10, 6, 2 
        expected_weight_10 =  [(0, 1)]
        expected_weight_6 =  [(0, 3), (1, 2)]
        expected_weight_2 =  [(0, 5), (1, 4), (2, 3)]
        self.assertEqual(get_fixed_weight(n,m,10,space),expected_weight_10)

        self.assertEqual(get_fixed_weight(n,m,6,space),expected_weight_6)

        self.assertEqual(get_fixed_weight(n,m,2,space),expected_weight_2)

        expected_weight_10 = [(0, 1)]
        self.assertEqual(get_fixed_weight(n,m,10,space),expected_weight_10)
        expected_weight_6 = [(0, 3), (1, 2)]
        self.assertEqual(get_fixed_weight(n,m,6,space),expected_weight_6)
        expected_weight_2 = [(0, 5), (1, 4), (2, 3)]
        self.assertEqual(get_fixed_weight(n,m,2,space),expected_weight_2)

        expected_weight_12= []
        self.assertEqual(get_fixed_weight(n,m,12,space),expected_weight_12)
        expected_weight_8= [(0, 2)]
        self.assertEqual(get_fixed_weight(n,m,8,space),expected_weight_8)
        expected_weight_4= [(0, 4), (1, 3)]
        self.assertEqual(get_fixed_weight(n,m,4,space),expected_weight_4)

        # in elements of weight 10
        self.assertEqual(rho_E(n,m,space,(0, 1)), [] )
        # in elements of weight 6
        self.assertEqual(rho_E(n,m,space,(0, 3)), [[12, 0, 2]] )
        self.assertEqual(rho_E(n,m,space,(1, 2)), [[6, 0, 2]] )
        # in elements of weight 2
        self.assertEqual(rho_E(n,m,space,(0, 5)), [[10, 0, 4]] )
        self.assertEqual(rho_E(n,m,space,(1, 4)), [[6, 0, 4], [12, 1, 3]] )
        self.assertEqual(rho_E(n,m,space,(2, 3)), [[10, 1, 3]] )
    

        # Check kernel of rho_E
        expected_ker_w_10 =  [1]
        expected_ker_w_6 =  [ 1, -2]
        expected_ker_w_2 =  [ 3, -5,  6]

        self.assertListEqual(kernel_rho_E(expected_weight_10,expected_weight_12,n,m,space),expected_ker_w_10)
        self.assertListEqual(kernel_rho_E(expected_weight_6,expected_weight_8,n,m,space),expected_ker_w_6)
        self.assertListEqual(kernel_rho_E(expected_weight_2,expected_weight_4,n,m,space),expected_ker_w_2)

        # Check rho_F
        # weight  10
        expected_rho_F_01 =  [[1, 0, 2]]
        self.assertListEqual(rho_F(n,m,space,(0, 1)),expected_rho_F_01)
        # weight  6
        expected_rho_F_03 =  [[1, 1, 3], [1, 0, 4]]
        self.assertListEqual(rho_F(n,m,space,(0, 3)),expected_rho_F_03)
        expected_rho_F_12 =  [[1, 1, 3]]
        self.assertListEqual(rho_F(n,m,space,(1, 2)),expected_rho_F_12)
        # weight  2
        expected_rho_F_05 =  [[1, 1, 5], [1, 0, 6]]
        self.assertListEqual(rho_F(n,m,space,(0, 5)),expected_rho_F_05)
        expected_rho_F_14 =  [[1, 2, 4], [1, 1, 5]]
        self.assertListEqual(rho_F(n,m,space,(1, 4)),expected_rho_F_14)
        expected_rho_F_23 =  [[1, 2, 4]]
        self.assertListEqual(rho_F(n,m,space,(2, 3)),expected_rho_F_23)

        # Check linear_rho_F
        # weight  10
        u_10 =  [[1, 0, 1]]
        expected_rho_w_10 =  [[1, 0, 2]]
        self.assertListEqual(linear_rho_F(u_10,n,m,space),expected_rho_w_10)
        # weight  6
        u_6 =  [[1, 0, 3], [-2, 1, 2]]
        expected_rho_w_6 =  [[-1, 1, 3], [1, 0, 4]]
        self.assertListEqual(linear_rho_F(u_6,n,m,space),expected_rho_w_6)
        # weight  2
        u_2 =  [[3, 0, 5], [-5, 1, 4], [6, 2, 3]]
        expected_rho_w_2 =  [[-2, 1, 5], [3, 0, 6], [1, 2, 4]]
        self.assertListEqual(linear_rho_F(u_2,n,m,space),expected_rho_w_2)

        # Check rho_F_iterates
        expected_iterates_w_10 =  [[[1, 0, 1]], [[1, 0, 2]], [[1, 1, 2], [1, 0, 3]], [[2, 1, 3], [1, 0, 4]], [[2, 2, 3], [3, 1, 4], [1, 0, 5]], [[5, 2, 4], [4, 1, 5], [1, 0, 6]], [[5, 3, 4], [9, 2, 5], [5, 1, 6]], [[14, 3, 5], [14, 2, 6]], [[14, 4, 5], [28, 3, 6]], [[42, 4, 6]], [[42, 5, 6]]]
        self.assertListEqual(rho_F_iterates(u_10,n,m,space,10),expected_iterates_w_10)
        expected_iterates_w_6 =  [[[1, 0, 3], [-2, 1, 2]], [[-1, 1, 3], [1, 0, 4]], [[-1, 2, 3], [1, 0, 5]], [[-1, 2, 4], [1, 1, 5], [1, 0, 6]], [[-1, 3, 4], [2, 1, 6]], [[-1, 3, 5], [2, 2, 6]], [[-1, 4, 5], [1, 3, 6]]]
        self.assertListEqual(rho_F_iterates(u_6,n,m,space,6),expected_iterates_w_6)
        expected_iterates_w_2 =  [[[3, 0, 5], [-5, 1, 4], [6, 2, 3]], [[-2, 1, 5], [3, 0, 6], [1, 2, 4]], [[-1, 2, 5], [1, 1, 6], [1, 3, 4]]]
        self.assertListEqual(rho_F_iterates(u_2,n,m,space,2),expected_iterates_w_2)

        # Check gen_repr
        expected_repr_10 = [[[1, 0, 1]], [[1, 0, 2]], [[1, 1, 2], [1, 0, 3]], [[2, 1, 3], [1, 0, 4]], [[2, 2, 3], [3, 1, 4], [1, 0, 5]], [[5, 2, 4], [4, 1, 5], [1, 0, 6]], [[5, 3, 4], [9, 2, 5], [5, 1, 6]], [[14, 3, 5], [14, 2, 6]], [[14, 4, 5], [28, 3, 6]], [[42, 4, 6]], [[42, 5, 6]]]
        self.assertListEqual(gen_repr(n,m,10,space), expected_repr_10)
        expected_repr_6 = [[[1, 0, 3], [-2, 1, 2]], [[-1, 1, 3], [1, 0, 4]], [[-1, 2, 3], [1, 0, 5]], [[-1, 2, 4], [1, 1, 5], [1, 0, 6]], [[-1, 3, 4], [2, 1, 6]], [[-1, 3, 5], [2, 2, 6]], [[-1, 4, 5], [1, 3, 6]]]
        self.assertListEqual(gen_repr(n,m,6,space), expected_repr_6)
        expected_repr_2 = [[[3, 0, 5], [-5, 1, 4], [6, 2, 3]], [[-2, 1, 5], [3, 0, 6], [1, 2, 4]], [[-1, 2, 5], [1, 1, 6], [1, 3, 4]]]
        self.assertListEqual(gen_repr(n,m,2,space), expected_repr_2)

if __name__=='__main__':
    from submodule_generation import *
    unittest.main()