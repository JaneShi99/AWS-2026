load("/home/janeshi/AWS26/Harvey-Class/P7-prec-3/P7-utils.sage")


def construct_T_bar_ijp(i, j, d, F_coeffs):
    arr = [[ZP3(0, 0, 0) for _ in range(d)] for _ in range(d)]

    for k in range(1, d):
        arr[k][k-1] = ZP3(2*j*F_coeffs[0], 2*i*F_coeffs[0], 0)
    
    for k in range(d):
        arr[k][d-1] = ZP3((d-k-2*j)*F_coeffs[d-k], (d-k-2*i)*F_coeffs[d-k], 0)

    return ZP3Matrix.from_array(ZZ, arr)


'''
Compute the single matrix T_k
'''
def generate_T_matrix(p, k_val, d, m, F_coeffs, R):
    sub = R(k_val) * F_coeffs[0]
    M = matrix(R, d, d)
    for i in range(1,d):
        M[i,i-1] = sub
    
    for i in range(d):
        j = d - i
        M[i, d-1] = (R(j*(m+1)) - R(k_val)) * F_coeffs[j]

    return M


# custom division in Z/p^mu Z
def divide_custom(x, y, p, R):
    divide_p_power = valuation(Integer(y), p)
    x_div_p = R(Integer(x)/(p^divide_p_power))
    y_div_p = R(Integer(y)/(p^divide_p_power))
    ans = R(x_div_p*(y_div_p^(-1)))
    return ans


def compute_A_f_avg_poly(F_coeffs, N):

    d = len(F_coeffs) - 1
    g = (d-1) // 2

    p_to_mat = [{}, {}]

    p_to_int = [{}, {}]


    # make the acc. remainder tree for matrices
    for i in [0, 1]:
        value_tree_leaves = []

        for j in range(1, N+1):
            T_bar_ijp = construct_T_bar_ijp(i, j, d, F_coeffs)
            value_tree_leaves.append(T_bar_ijp)

        value_tree = build_product_tree(value_tree_leaves)
        modulus_tree = build_product_tree([k^2 if is_prime(k) else 1 for k in range(1, N+1)])

        rem_tree, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZP3Matrix.identity(ZZ, d))
        
        '''
        for x, y in enumerate(leaf_val_list):
            print(x, y)
        '''
        
        for p in range(N+1):
            if is_prime(p) and p != 2:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                T_bar_mod_p3 = T_bar_int.apply_map(lambda x: Zmod(p^3)(x))
                T_mod_p3 = T_bar_mod_p3.apply_map(lambda x: Zmod(p^3)(x * inverse_mod(2^(p-1), p^3)))
                p_to_mat[i][p] = T_mod_p3
    
    #make the acc remainder tree for integers
    for i in [0, 1]:
        value_tree_leaves = []

        for k in range(1, N+1):
            int_i_j = ZP3(k, i, 0)
            value_tree_leaves.append(int_i_j)
        
        value_tree = build_product_tree(value_tree_leaves)
        modulus_tree = build_product_tree([k^3 if is_prime(k) else 1 for k in range(1, N+1)])

        rem_tree, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZP3(1,0,0))

        for p in range(N+1):
            if is_prime(p) and p != 2:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                p_to_int[i][p] = Zmod(p^3)(T_bar_int)
    
    

    #now, we're ready to compute A_f fully. We have a loop that computes A_f for each p at a time

    p_to_A_f = {}

    for p in range(N+1):
        if is_prime(p) and p > 5:

            R = Zmod(p^3)
            m = (p-1)//2

            F_coeffs_p3 = [R(c) for c in F_coeffs]
            F0_p3 = R(F_coeffs_p3[0])

            acc = Matrix(R, 1, d)

            R_0 = p_to_mat[0][p]
            R_p = p_to_mat[1][p]

            Rp_minus_R0 = R_p - R_0
            R_1 = Rp_minus_R0.apply_map(lambda x: R(Integer(x)/p))
            sprint_matrices = [R_0 + (l*p)*R_1 for l in range(0, g)]

            int_0 = p_to_int[0][p]
            int_p = p_to_int[1][p]

            int_p_minus_int_0 = int_p - int_0
            int_1 = R(Integer(int_p_minus_int_0)/p)

            int_products = [int_0 + (l*p)*int_1 for l in range(0, g)]

            A_f = []

            for l in range(0, g):
                # the following step is to compute the single step 
                # U_(ip) = U_(ip-1)*T_(ip)*(constants)
                if l == 0:
                    acc[0, -1] = R(F0_p3^m)
                else:
                    T_single = generate_T_matrix(p, p*l, d, m, F_coeffs, R)

                    Hk = acc * T_single[:,-1]
                    to_invert = l*p*F0_p2 
                    Hk = Hk.apply_map(lambda x: divide_custom(x, to_invert, p, R))
                    acc_new = [acc[0][i] for i in range(1, len(acc[0]))] + [Hk[0][0]]
                    acc = Matrix(R, 1, d, acc_new)

                acc = acc * sprint_matrices[l]
                to_invert = int_products[l]*(F0_p3^(p-1))
                acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, R))
                
                A_f.append(list(acc[0][(-g):]))

            A_f_p = Matrix(A_f).apply_map(lambda x: Integer(x) % p)
            p_to_A_f[p] = A_f_p
    
    return p_to_A_f



'''
AVG poly up to 50,000 took 59 seconds
Sqrt up to 50,000 took 1943 seconds (30 minutes) 
'''
        
    
start = timer()
#N = 102
N = 6300
R.<x> = PolynomialRing(Integers())
#f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
f =  -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
answer = compute_A_f_avg_poly(f.list(), N)
for key, value in answer.items():
    print(key)
    print(value)


time = timer() - start
print(time)