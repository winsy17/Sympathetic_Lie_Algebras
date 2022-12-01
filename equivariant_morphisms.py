from fractions import Fraction
import warnings
import copy
import math
from submodule_generation import gen_repr

def get_basis(n,m,space):
    '''
    Input: integers n and m such that 1 <= n <= m,  
    Output: It returns a list with elements of a basis of 
    V_n tensor V_m. The elements of this list are tuples (i,j)
    corresponding to elements e_i and f_j of the basis of V_n and V_m respectively
    '''
    if n > m:
        raise ValueError('Arguments n, m should satisfy n <= m')
    extra = 0
    c = 0
    if space != 'tensor':
        c = 1
        if space == 'exterior power':
            extra = 1
    return([(i,j) for i in range(n+1) for j in range(c*i+extra,m+1)])

def get_limits(n,space):
    low_lim = 0
    upp_lim = math.floor(n/2)
    if space == 'exterior power':
        low_lim += 1
        upp_lim = math.floor((n+1)/2)
    elif space == 'tensor':
        upp_lim = n
    return(low_lim,upp_lim)

def get_index(q,space):
    if space == 'tensor':
        return(q)
    elif space == 'exterior power':
        return(2*q-1) 
    elif space == 'symmetric power':
        return(2*q)

def det(M):
    '''
    This function was written by 'pycoder' and extracted from: 
    https://stackoverflow.com/questions/66192894/precise-determinant-of-integer-nxn-matrix
    '''
    M = [row[:] for row in M] # make a copy to keep original M unmodified
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0: # swap with another row having nonzero i's elem
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return(0) # all M[*][i] are zero => zero determinant
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                assert ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) % prev == 0
                M[j][k] = ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) // prev
        prev = M[i][i]
    return(sign * M[-1][-1]) 

def det_minors(Mat,i,j):
    ''' Input: M is a list of lists representing a square matrix with integer entries
    i,j are integers 0<= i,j <= len(M)
    output: the determinant of the (i,j)-th first minor 
    of a matrix M
    '''
    minor = copy.deepcopy(Mat)
    del minor[i]
    return(det([row[:j] + row[j+1:] for row in minor]))


def get_inverse(Mat):
    '''Input: Mat is a list representing a squared matrix M of integer entries
    Output: list Inv_Mat, the transpose of the matrix of cofactors of M'''
    inv_det = Fraction(1,det(Mat))
    if len(Mat) == 1:
        return([[inv_det]])
    Inv_Mat= [[inv_det*((-1)**(i+j))*det_minors(Mat,i,j) for i in range(len(Mat))] for j in range(len(Mat))]
    return(Inv_Mat)


def equiv_morphism(n, m, k):
    '''
    Input: integers n and m such that 1 <= n <= m,  
    postive integer k
    Output: It returns a list of dimensions ((n+1)*(m+1),(k+1))
    This array is the matrix associated to the equivariant morphism
    in terms a basis B of V_n tensor V_m (calculated below)
    If n+m-k mod 2 is non zero, it returns a zero matrix of the same dimensions. 
    '''
    
    # Raise input error if n > m
    if n > m:
        message = 'Expected n <= m'
        raise ValueError(message)
    
    # Check n+m-k=0 mod2 and m-n <= k <= n+m
    if (n+m-k)%2!= 0 or k < (m-n) or (n+m) < k:
        message = 'The only bilinear sl_2(C)-equivariant map T:V_'+str(n)+' x V_'+str(m)+' --> V_'+str(k)+ ' is trivial'
        warnings.warn(message)
        return([[0 for i in range(k+1)] for j in range((n+1)*(m+1))])

    # When n==m, the value of k completely determines whether the morphism is symmetric or skew-symmetric    
    if n == m: 
        if (2*n-k)%4 == 0:
            space = 'symmetric power'
        else:# (2*n-k)%4 == 2:
            space = 'exterior power'
    else:
        # else construct a bilinear morphism
        space = 'tensor'
    
    ##Obtain basis
    B = get_basis(n,m,space)
    Mat = []
    low_lim, upp_lim = get_limits(n,space)
    index_k = 0

    ##Get basis of V_k, for each term V_k in the Clebsh-Gordan decomposition 
    for q in range(low_lim, upp_lim+1):
        p = n+m-2*get_index(q,space)
        V_p = gen_repr(n,m,p,space)
        # writing the basis in V_p in terms of the basis of tensor prod/ext. power/sym. power
        for v in V_p:
            v_B = [0 for i in range(len(B))]
            for [c,i,j] in v:
                if space in ['symmetric power','exterior power'] and i != j:
                    v_B[B.index((i,j))] = 2*int(c)
                else:
                    v_B[B.index((i,j))] = int(c)
            Mat.append(v_B)
        if p > k:
            index_k += p + 1 

    ## Getting inverse of Mat by inverting blocks of fixed weight.
    inv = [[0 for i in range(len(Mat))] for j in range(len(Mat))]
    if space == 'exterior power':
        c = 1
    else:
        c = 0
    ## Looping through all possible weights
    for s in range(c,n+m-c+1):
        col = [B.index(b) for b in B if b[0] + b[1] == s]
        s_aux = s 
        if (n+m-2*s) < 0:
            s_aux = n+m-s
        
        # Looping through all the V_k's to find basis elements
        # with weight = n+m-2*s.
        # V_k contains a basis element with weight n+m-2*s
        # if and only if k >= |n+m-2*s|
        if space == 'tensor':
            row = [s + a*(n+m) + a*(1-a) for a in range(min(s_aux,n)+1)]
        elif space == 'symmetric power':
            row = [s + 2*n*a + a*(1-2*a) for a in range(math.floor(s_aux/2)+1)]
        else:
            row = [s + 2*(a-1)*(n+1-a) - a for a in range(1,math.floor((s_aux+1)/2)+1)]
        
        Mat_w = [[Mat[r][c] for c in col] for r in row]
        Inv_w = get_inverse(Mat_w)

        for r in range(len(row)):
            for c in range(len(col)):
                inv[col[c]][row[r]] = Inv_w[c][r]

    #Obtaining a solution for the linear system
    temp = [row[index_k:index_k+k+1] for row in inv]
    Lcm = 1
    for i in range(len(temp)):
        for j in range(len(temp[0])):
            if temp[i][j] !=0:
                Lcm = math.lcm(Lcm,temp[i][j].denominator)
    prod = [[0 for i in range(len(temp[0]))] for j in range(len(temp))]
    for i in range(len(temp)):
        for j in range(len(temp[0])):
            prod[i][j] = int(Lcm*temp[i][j])    
    if space == 'tensor':
        return(prod)
    else:
        # Obtain basis of tensor product
        Basis = get_basis(n,m,'tensor')
        # Extending a morphism into a symmetric/skew-symmetric morphism on the tensor product
        morphism = [[0 for i in range(k+1)] for j in range(len(Basis))]
        sgn = 1
        if space == 'exterior power':
            sgn = -1
        for (i,j) in B:
            morphism[Basis.index((i,j))] = [f for f in prod[B.index((i,j))]]  
            morphism[Basis.index((j,i))] = [sgn*f for f in prod[B.index((i,j))]]      
        return(morphism)

if __name__ == "__main__":
    n = 6
    m = 6
    k = 2

    A = equiv_morphism(n, m, k)