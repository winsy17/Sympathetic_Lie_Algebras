import itertools
from sympy import *

def left_mult(A): #DA
    '''
    A is a list of lists, representing a matrix of nxn n=len(A)
    Output: returns a list containing the equations of DA (where D is the unknown matrix of nxn - vector of n^2)
    '''
    L_mult = []
    n= len(A)
    for i in range(n):
        for column in zip(*A):
            mult =  [0 for r in range(i*(n))]+[elem for elem in column] + [0 for r in range((n-i-1)*n)]
            L_mult.append(mult)
    return(L_mult)

def right_mult(A): # (-1)*AD
    '''
    A is a list of lists, representing a matrix of nxn n=len(A)
    Output: returns a list containing the equations of AD (where D is the unknown matrix of nxn - vector of n^2)
    '''
    R_mult = []
    n= len(A)
    c = -1
    for row in A:
        for i in range(n):
            mult = []
            for j in range(n):
                mult += [0 for r in range(i)] + [c*row[j]] + [0 for r in range(n-i-1)]
            R_mult.append(mult)
    return(R_mult)

def mix_adj(ind,ad): # (-1)*(sum over k of d_{ki}ad_(e_k))
    ''' ind integer 0 <= ind < n, where n is the dimension of the Lie algebra
        ad is a list of n adjoint matrices of n by n
    '''
    c= -1
    n = len(ad)
    L_mix = []
    for i in range(n):
        for j in range(n):
            v = [0 for r in range(n*n)]
            for k in range(n):
                v[ind + k*n] = c*ad[k][i][j]
            L_mix.append(v)
    return(L_mix)

def vector_sum(a,b):
    '''
    a and b are lists of the same length
    output: the term-wise sum of these lists
    '''
    return([a[i] + b[i] for i in range(len(a))])

def sum_terms(A,B):
    '''
    A,B lists of lists each sublist has length n^2
    '''
    sum_lists = []
    for i in range(len(A)):
        sum_lists.append(vector_sum(A[i],B[i]))
    return(sum_lists)

def scalar_mult(c,A):
    '''c a number
    A a list of lists
    Each sublist gets multiplied by c'''
    cA = []
    for row in A:
        mult_c = [elem*c for elem in row]
        cA.append(mult_c)
    return(cA)

def space_of_derivations(ad):
    '''
    ad is a list of n adjoint matrices 
    output: basis for the vector space of derivations
    '''
    n = len(ad)
    Mat = []
    for i in range(n):
        # Get equations for D*ad_{e_i} - sum_k d_{ki}*ad_{e_k} - ad_{e_i}*D
        Mat += sum_terms(sum_terms(left_mult(ad[i]),mix_adj(i,ad)),right_mult(ad[i]))
    
    Mat.sort()
    Mat_simple = list(Mat for Mat,_ in itertools.groupby(Mat))

    flat_ad = Matrix(1,n*n,Matrix(ad[0]))
    
    for i in range(1,len(ad)):
        flat_ad = flat_ad.row_insert(i+1,Matrix(1,n*n,Matrix(ad[i])))
    
    ker = Matrix(Mat_simple).nullspace()
    ker_t = [Matrix(n,n,k.transpose()) for k in ker]

    print('The space of derivations has dimension ' + str(len(ker_t)))
    return(ker_t)


if __name__ == "__main__":
    ad_H = [[0, 0, 0], [0, 2, 0], [0, 0, -2]]
    ad_E = [[0, 0, 1], [-2, 0, 0], [0, 0, 0]]
    ad_F = [[0, -1, 0], [0, 0, 0], [2, 0, 0]]
    adjoints = [ad_H, ad_E, ad_F]
    ker = space_of_derivations(adjoints)
    print(ker)