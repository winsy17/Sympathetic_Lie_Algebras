import math
from fractions import Fraction


def get_fixed_weight(n,m,k,space):
    '''Input:
    integers - n, m, k, n<=m
    string - space: 'tensor', 'symmetric power', 'exterior power'.
    Output: A list of indices (i,j) such that (e_i,f_j)
    is an element of weight k. This is calculated using the representation rho(H)'''
    if (n+m-k)%2 != 0:
        raise ValueError('There are no elements with weight '+ str(k) + ' in the tensor product of V_' +str(n)+ ' with V_' +str(m))
    if space in ['symmetric power', 'exterior power'] and n != m:
        raise ValueError('Expected n = m for exterior/symmetric power')
    extra = 0
    c = 0
    a = int((n+m-k)/2)
    if space != 'tensor':
        c = 1
        if space == 'exterior power':
            extra = 1
    return([(i,a-i) for i in range(n+1) if (c*i+extra)<= a-i <= m])

def rho_E(n,m,space,t):
    '''Input:
        Integers n <= m to calculate rho_n tensor rho_m (E)
        space is a string : 'tensor', 'symmetric power', 'exterior power'.
        The last two inputs of space only make sense when n=m.
        t  is a tuple of the form (i,j) containing the indices of the elements e_i tensor f_j
        Output: A list with elements of the form [c,r,s], which corresponds
        to the factor of rho_E of the form c(e_r tensor f_s). 
        If rho_E(e_i tensor f_j) = 0, it returns an empty list.
    '''
    i = t[0]
    j = t[1]
    if space == 'symmetric power' and i>j:
        raise ValueError("Expected i<= j for symmetric power") 
    if space == 'exterior power' and i>=j:
        raise ValueError("Expected i<j for exterior power") 
    if i == j and i > 0 and space == 'symmetric power':
        return([[i*(n+1-i),i-1,i]]) # if e_i tensor e_i is the element of the basis
        # return([[2*i*(n+1-i),i-1,i]]) # if 2*(e_i tensor e_i) is the element of the basis
    if i == j-1 and space == 'exterior power':
        if i > 0:
            return([[i*(n+1-i),i-1,i+1]])
        else:
            return([])
    image = []
    if i > 0:
        image.append([i*(n+1-i),i-1,j])
    if j > 0:
        image.append([j*(m+1-j),i,j-1])
    return(image)

def kernel_rho_E(W_k,W_k2,n,m,space):
    '''
    Input:
    W_k: list of tuples, a tuple (i,j) represents 
    the tensor of e_i and f_j, which is an element of weight k
    W_k2: list of tuples of weight k+2
    n <= m integers
    string space: 'tensor', 'symmetric power', 'exterior power'.
    Output: returns the element of degree k in the kernel of rho_E
    as a linear combination of the elements of W_k2 (in the order listed in W_k2)
    '''
    if W_k2 == []:
        return([1 for i in range(len(W_k))])
    matrix = []
    for t in W_k:
        rho = rho_E(n,m,space,t)
        mat = [0 for i in range(len(W_k2))]#np.zeros(len(W_k2))
        for (c,r,s) in rho:
            l = W_k2.index((r,s))
            mat[l] = c
        matrix.append(mat)
    coeff = [Fraction(1,1)]
    Lcm = 1
    for i in range(1,len(W_k)):
        coeff+=[Fraction((-1)*int(matrix[i-1][i-1]),int(matrix[i][i-1]))*coeff[i-1]]
        Lcm = math.lcm(Lcm,coeff[i].denominator)
    ker_generator = [int(Lcm*c) for c in coeff]
    return(ker_generator)

def rho_F(n,m,space,t):
    '''Input:
    Integers n <= m to calculate rho_n tensor rho_m (F)
        space is a string : 'tensor', 'symmetric power', 'exterior power'.
        The last two options for space only make sense when n=m.
        t  is a tuple of the form (i,j) containing the indices of the elements e_i tensor f_j
        Output: A list with elements of the form [c,r,s], which corresponds
        to the factor of rho_F of the form c(e_r tensor f_s). 
        If rho_F(e_i tensor f_j) = 0, it returns an empty list.
    '''
    i = t[0]
    j = t[1]
    if space == 'symmetric power' and i>j:
        raise TypeError("Error: expected i<= j for symmetric power") 
    if space == 'exterior power' and i>=j:
        raise TypeError("Error: expected i<j for exterior power") 
    image = []
    if space == 'symmetric power' and i == j and i < n:
        return([[1,i,i+1]]) # e_i tensor e_{i+1} + e_{i+1} tensor e_i 
    if i < n:
        if i != j-1:
            image.append([1,i+1,j])
        elif space == 'tensor':
            image.append([1,i+1,j])
        elif space == 'symmetric power':
            image.append([2,i+1,j])
    if j < m:
        image.append([1,i,j+1])
    return(image) 

def linear_rho_F(u,n,m,space):
    '''Input:
        A list u of triples [c,r,s] representing
        a linear combination of elements c*(e_r tensor f_s), c is a scalar.
        Integers n <= m to calculate rho_n tensor rho_m (F)
        space is a string : 'tensor', 'symmetric power', 'exterior power'.
        The last two options for space only make sense when n=m.
        Output: This function evaluates the linear map rho_F on u. 
        It returns a list containing the summands of rho_F(u)
        If rho_F(e_i tensor f_j) = 0, it returns an empty list.
    '''
    F = []
    for [c,i,j] in u:
        R = rho_F(n,m,space,(i,j)) # this is a list
        for [x,y,z] in R:
            indicator = 0
            for r in range(len(F)):
                if y == F[r][1] and z == F[r][2]:
                    F[r][0] += c*x
                    indicator = 1
            if indicator == 0:
                F.append([c*x,y,z])
    simple_F = [[c,i,j] for [c,i,j] in F if c!=0]
    return(simple_F)

def rho_F_iterates(u,n,m,space,k):
    '''Input:
        A list u of triples [c,r,s] which corrensponds to the 
        element c*(e_r tensor f_s), c is a scalar.
        Integers n <= m to calculate rho_n tensor rho_m (F)
        space is a string : 'tensor', 'symmetric power', 'exterior power'.
        The last two options for space only make sense when n=m.
        k >= 1 is the number of iterates to calculate
        Output: Returns a list with [u,rho_n_m(F)(u),...,(rho_n_m(F)(u))^k]
    '''
    iterated = [u]
    for i in range(k):
        v = linear_rho_F(iterated[i],n,m,space)
        if v == []:
            raise TypeError("Error: input vector u is not generating the full subspace V_k")
        else:
            iterated.append(v)
    return(iterated)

def gen_repr(n,m,k, space = 'tensor'):
    '''
    Input: integers n and m such that 1 <= n <= m,  
           integer n+m-k is an even number
           string space: 'tensor', 'symmetric power', 'exterior power'.
           The last two inputs of space only make sense when n=m.
    Output: It returns a basis of the irrep V_k as a simple module of
    V_n tensor V_m or of the symmetric/antisymmetric power of V_n
    It returns an error if V_k is not a submodule of the decomposition
    '''
    # Raise input error if n > m
    if n > m:
        message = 'Expected n <= m'
        raise ValueError(message)
    if space != 'tensor' and n!= m:
        raise ValueError('Expected n=m for ' + space)

    # Raise error if V_k is not a submodule of the decomposition
    # Check n+m-k=0 mod2 and m-n <= k <= n+m
    if (n+m-k)%2!= 0 or k < (m-n) or (n+m) < k:
        if space == 'tensor':
            message = 'V_'+str(k)+' is not a submodule of the tensor of V_' + str(n) + ' with V_' + str(m)   
        elif space == 'symmetric power':
            message = 'V_'+str(k)+' is not a submodule of S^2(V_' + str(n) + ')'
        else: 
            message = 'V_'+str(k)+' is not a submodule of Ext^2(V_' + str(n) + ')'
        raise ValueError(message)

    # Check (2*n-k) mod 4
    if space == 'symmetric power' and (2*n-k)%4 != 0:
        raise ValueError('V_'+str(k)+' is not a submodule of S^2(V_' + str(n) + ')')
    if space == 'exterior power' and (2*n-k)%4 != 2:
        raise ValueError('V_'+str(k)+' is not a submodule of Ext^2(V_'+ str(n) + ')')
    
    # Obtain a non-zero element of weight k
    # in the kernel of the tensor representation

    # Obtain vectors of weight k
    W_k = get_fixed_weight(n,m,k,space)

    # Obtain vectors of weight k+2 - rho_E increases weight by 2
    W_k2 = get_fixed_weight(n,m,k+2,space)
    
    # Obtain non zero vector of degree k, v, of kernel of rho_E
    # v will be expressed as a linear combination of elements in W_k
    # WARNING u will be in the same order than W_k rather than the order of the basis L
    if W_k2 == []:
        u = []
        u = [[1,W_k[i][0],W_k[i][1]] for i in range(len(W_k))] 
    else:
        v = kernel_rho_E(W_k,W_k2,n,m,space)
        # Transform v into the input format that rho_F uses:
        # [c_0,...,c_r] to [[c_0,i_0,j_0],...,[c_r,i_r,j_r]]
        u = []
        u = [[v[i],W_k[i][0],W_k[i][1]] for i in range(len(W_k))] 

    # Obtain the iterates of rho_F
    span = rho_F_iterates(u,n,m,space,k)
    return(span)

if __name__ == "__main__":
    n = 3
    m = 3
    k = 6
    
    space = 'symmetric power'

    Rep= gen_repr(n,m,k,space)
    print(Rep)


    
    