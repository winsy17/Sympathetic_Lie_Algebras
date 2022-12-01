import numpy as np

def get_adjoints(struc_const):
    '''
    struc_const: list of length n(n-1)/2 that contains the 
    structure constants of the Lie algebra as lists 
    of length n.
    struc_const[k]=[c_0,...,c_n] corresponds a Lie bracket [e_i,e_j] = [c_0,...,c_{n-1}] = sum c_re_r
    with i<j = 0,...,n-1
    This function returns a list with the n adjoint matrices ad(e_i)
    '''
    n = len(struc_const[0])
    ad = []

    for k in range(n):
        E_k = [[(-1)*elem for elem in struc_const[j*n-(j*(j+1))//2+k-j-1]] for j in range(k)] + [[0 for r in range(n)]] +[struc_const[n*k-(k*(k+1))//2+j-k-1] for j in range(k+1,n)]
        ad+=[list(map(list, zip(*E_k)))]
    return(ad)


def is_Lie_algebra(struc_const):
    '''
    struc_const: list of length n(n-1)/2 that contains the 
    structure constants of the Lie algebra as lists 
    of length n.
    struc_const[k]=[c_0,...,c_n] corresponds a Lie bracket [e_i,e_j] = [c_0,...,c_{n-1}] = sum c_re_r
    with i<j = 0,...,n-1
    This function checks if the Jacobi identity is satisfied and 
    returns True if the vector space <e_0,...,e_{n-1}> is a Lie with the above structure constants.
    '''
    n = len(struc_const[0])
    ad = get_adjoints(struc_const)
    e = np.identity(n,dtype = int)
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                ad_ijk = np.matmul(np.array(ad[i]),np.matmul(np.array(ad[j]),e[k]))
                ad_jki = np.matmul(np.array(ad[j]),np.matmul(np.array(ad[k]),e[i]))
                ad_kij = np.matmul(np.array(ad[k]),np.matmul(np.array(ad[i]),e[j]))
                J_ijk = ad_ijk + ad_jki + ad_kij
                if np.amax(np.absolute(J_ijk)) != 0:
                    return(bool(False))    
    return(bool(True))
    


if __name__ == "__main__":
    str_const = [[0,2,0],[0,0,-2],[1,0,0]]
    print(is_Lie_algebra(str_const))