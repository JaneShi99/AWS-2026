# We wish to compute (p-1)! mod p^2

# calculating remainder tree from bottom up.
# scrap this! can be generalized

'''
def compute_Q_k(p,t):
    #compute (k+1)(k+2)...(k+t) mod p^2
    R.<k> = PolynomialRing(Zmod(p^2))

    def recursive_helper(a,b):
        if a == b:
            return R(1)
        elif a == b-1:
            return R(k+a)
        else:
            m = a + (b-a)//2
            return recursive_helper(a,m) * recursive_helper(m,b)
    
    return recursive_helper(1,t+1)
'''


def linear_poly_acc(R, lin_list):
    def recursive_helper(a,b):
        if a == b:
            return R(1)
        elif a == b-1:
            return lin_list[a]
        else:
            m = a + (b-a)//2
            return recursive_helper(a,m) * recursive_helper(m,b)
    
    return recursive_helper(0, len(lin_list))

def compute_Q_k(p, t):
    R.<k> = PolynomialRing(Zmod(p^2))
    lin_list = [k+i for i in range(1, t+1)]
    return linear_poly_acc(R, lin_list)

print(compute_Q_k(11, 16))


#calculating remainder tree from top down.

