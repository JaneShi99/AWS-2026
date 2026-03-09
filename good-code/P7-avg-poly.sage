import os as _os
if _os.path.exists("P7-utils.sage"):
    load(_os.path.abspath("P7-utils.sage"))
else:
    load(_os.path.abspath("good-code/P7-utils.sage"))

# Constructing the matrix on page 39 of notes (Section 7.3)
def construct_T_bar_ijp(i, j, d, F_coeffs):
    arr = [[ZP2(0, 0) for _ in range(d)] for _ in range(d)]

    for k in range(1, d):
        arr[k][k-1] = ZP2(2*j*F_coeffs[0], 2*i*F_coeffs[0])
    
    for k in range(d):
        arr[k][d-1] = ZP2((d-k-2*j)*F_coeffs[d-k], (d-k-2*i)*F_coeffs[d-k])

    return ZP2Matrix.from_array(ZZ, arr)


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


def compute_A_f_avg_poly(F_coeffs, d, N):
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

        rem_tree, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZP2Matrix.identity(ZZ, d))
        

        for p in range(N+1):
            if is_prime(p) and p != 2:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                T_bar_mod_p2 = T_bar_int.apply_map(lambda x: Zmod(p^2)(x))
                T_mod_p2 = T_bar_mod_p2.apply_map(lambda x: Zmod(p^2)(x * inverse_mod(2^(p-1), p^2)))
                p_to_mat[i][p] = T_mod_p2
    
    #make the acc remainder tree for integers
    for i in [0, 1]:
        value_tree_leaves = []

        for k in range(1, N+1):
            int_i_j = ZP2(k, i)
            value_tree_leaves.append(int_i_j)
        
        value_tree = build_product_tree(value_tree_leaves)
        modulus_tree = build_product_tree([k^2 if is_prime(k) else 1 for k in range(1, N+1)])

        rem_tree, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZP2(1,0))

        for p in range(N+1):
            if is_prime(p) and p != 2:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                p_to_int[i][p] = Zmod(p^2)(T_bar_int)
    
    

    #now, we're ready to compute A_f fully. We have a loop that computes A_f for each p at a time

    p_to_A_f = {}

    for p in range(N+1):
        if is_prime(p) and p > 5:
            # Skip bad primes: p divides F_coeffs[0] means f has a root at x=0 mod p
            if gcd(p, ZZ(F_coeffs[0])) != 1:
                continue
            # doesn't work if p < g
            if p < g:
                continue

            R = Zmod(p^2)
            m = (p-1)//2

            F_coeffs_p2 = [R(c) for c in F_coeffs]
            F0_p2 = R(F_coeffs_p2[0])

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
            
            '''
            for mat in sprint_matrices:
                print(mat)
                print("\n")
            print("int")
            print(int_products)
            '''

            for l in range(0, g):
                # the following step is to compute the single step 
                # U_(ip) = U_(ip-1)*T_(ip)*(constants)
                if l == 0:
                    acc[0, -1] = R(F0_p2^m)
                else:
                    T_single = generate_T_matrix(p, p*l, d, m, F_coeffs, R)

                    Hk = acc * T_single[:,-1]
                    to_invert = l*p*F0_p2 
                    Hk = Hk.apply_map(lambda x: divide_custom(x, to_invert, p, R))
                    acc_new = [acc[0][i] for i in range(1, len(acc[0]))] + [Hk[0][0]]
                    acc = Matrix(R, 1, d, acc_new)

                acc = acc * sprint_matrices[l]
                to_invert = int_products[l]*(F0_p2^(p-1))
                acc = acc.apply_map(lambda x: x*to_invert^(-1))
                # reversed: acc entries are ordered right-to-left in Harvey's
                # column convention, so the last entry corresponds to column 0
                A_f.append(list(reversed(acc[0][(-g):])))

            A_f_p = Matrix(A_f).apply_map(lambda x: Integer(x) % (p))
            p_to_A_f[p] = A_f_p
    
    return p_to_A_f


def compute_A_f_avg_poly_from_curve(C, N):
    F_coeffs_poly, _ = C.hyperelliptic_polynomials()
    d = C.degree()
    F_coeffs = [Integer(c.lift()) for c in F_coeffs_poly]
    return compute_A_f_avg_poly(F_coeffs, d, N)


'''
AVG poly up to 50,000 took 59 seconds
Sqrt up to 50,000 took 1943 seconds (30 minutes) 
'''
        
    
import os as _os
if 'P7-avg-poly' in _os.path.basename(sys.argv[0]):
    start = timer()
    N = 6300
    R.<x> = PolynomialRing(Integers())
    #f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
    f =  -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
    C = HyperellipticCurve(f)
    answer = compute_A_f_avg_poly_from_curve(C, N)
    for key, value in answer.items():
        print(key)
        print(value)

    time = timer() - start
    print(time)
