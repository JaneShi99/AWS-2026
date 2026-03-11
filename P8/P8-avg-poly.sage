import os as _os
from typing import reveal_type
load(_os.path.abspath("P8-utils.sage"))

# Constructing the matrix on page 39 of notes (Section 7.3)
def construct_T_bar_ijp(i, j, d, F_coeffs, mu, ell):
    arr = [[ZPmu([0, 0], mu) for _ in range(d)] for _ in range(d)]

    for k in range(1, d):
        arr[k][k-1] = ZPmu([2*j*F_coeffs[0], 2*i*F_coeffs[0]], mu)
    
    for k in range(d):
        arr[k][d-1] = ZPmu([(d-k-2*j)*F_coeffs[d-k], ((d-k)*(2*ell + 1) - 2*i)*F_coeffs[d-k]], mu)

    return ZPmuMatrix.from_array(ZZ, arr, mu)


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


def compute_A_F_l_avg_poly(F_coeffs, d, N, ell=0):
    g = (d - 1)//2
    lam: Integer = ceil(g/2)
    mu = lam + 1
    dl = (2*ell + 1)*(g + 1) - 1

    min_prime = max(4*g, ell*d + 1)

    p_to_mat = [{} for _ in range(dl)]
    p_to_int = [{} for _ in range(dl)]

    # make the acc. remainder tree for matrices
    for i in [0..(dl-1)]:
        value_tree_leaves = [ZPmuMatrix.identity(ZZ, d, mu) for _ in range(d*ell)]

        for j in [(d*ell + 1)..N]:
            T_bar_ijp = construct_T_bar_ijp(i, j - d*ell, d, F_coeffs, mu, ell)
            value_tree_leaves.append(T_bar_ijp)

        value_tree = build_product_tree(value_tree_leaves)
        modulus_tree = build_product_tree([k^mu if is_prime(k) else 1 for k in range(1, N+1)])

        _, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZPmuMatrix.identity(ZZ, d, mu))
        

        for p in range(N+1):
            if is_prime(p) and p > min_prime:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                T_bar_mod_p_mu = T_bar_int.apply_map(lambda x: Zmod(p^mu)(x))
                T_mod_p_mu = T_bar_mod_p_mu.apply_map(lambda x: Zmod(p^mu)(x * inverse_mod(2^(p-d*ell-1), p^mu)))
                p_to_mat[i][p] = T_mod_p_mu
    
    #make the acc remainder tree for integers
    for i in [0..(dl - 1)]:
        value_tree_leaves = [ZPmu([1], mu) for _ in range(d*ell)]

        for k in [(d*ell + 1)..N]:
            int_i_j = ZPmu([k, i], mu)
            value_tree_leaves.append(int_i_j)
        
        value_tree = build_product_tree(value_tree_leaves)
        modulus_tree = build_product_tree([k^mu if is_prime(k) else 1 for k in range(1, N+1)])

        _, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree, identity = ZPmu([1], mu))

        for p in range(N+1):
            if is_prime(p) and p > min_prime:
                T_bar = leaf_val_list[p-1]
                T_bar_int = T_bar.realize(p)
                p_to_int[i][p] = Zmod(p^mu)(T_bar_int)

    # compute the small-step matrices
    small_step_matrix = []
    for i in [0..(dl - 1)]: 
        i_th_row_small_step_matrices = []
        for s in [1..ell]:
            leaves = [
                construct_T_bar_ijp(i + 1, -s*d + j, d, F_coeffs, mu, ell)
                for j in range(d)
            ]
            i_th_row_small_step_matrices.append(prod(leaves))
        small_step_matrix.append(i_th_row_small_step_matrices)

    small_step_denominators = []
    for i in [0..(dl - 1)]: 
        i_th_row_small_step_denominators = []
        for s in [1..ell]:
            leaves = [
                ZPmu([-s*d + j, i + 1], mu)
                for j in range(d)
            ]
            i_th_row_small_step_denominators.append(prod(leaves))
        small_step_denominators.append(i_th_row_small_step_denominators)
        

    #now, we're ready to compute A_f fully. We have a loop that computes A_f for each p at a time

    p_to_A_f = {}

    for p in range(N+1):
        if is_prime(p) and p > min_prime:
            # Skip bad primes: p divides F_coeffs[0] means f has a root at x=0 mod p
            if gcd(p, ZZ(F_coeffs[0])) != 1:
                continue
            # doesn't work if p < g
            if p < g:
                continue

            R = Zmod(p^mu)
            m = (p-1)//2

            F_coeffs_p_mu = [R(c) for c in F_coeffs]
            F0_p_mu = R(F_coeffs_p_mu[0])

            acc = Matrix(R, 1, d)
            sprint_matrices = [p_to_mat[i][p] for i in range(0, dl)]
            int_products = [p_to_int[l][p] for l in range(0, dl)]

            A_f = []
            
            '''
            for mat in sprint_matrices:
                print(mat)
                print("\n")
            print("int")
            print(int_products)
            '''

            for l in range(0, dl):
                # the following step is to compute the single step 
                # U_(ip) = U_(ip-1)*T_(ip)*(constants)
                if l == 0:
                    acc[0, -1] = R(F0_p_mu^((2*ell + 1)*m))
                else:
                    T_single = generate_T_matrix(p, p*l, d, (2*ell + 1)*m, F_coeffs, R)

                    Hk = acc * T_single[:,-1]
                    to_invert = l*p*F0_p_mu 
                    Hk = Hk.apply_map(lambda x: divide_custom(x, to_invert, p, R))
                    acc_new = [acc[0][i] for i in range(1, len(acc[0]))] + [Hk[0][0]]
                    acc = Matrix(R, 1, d, acc_new)

                current_row = []
                acc = acc * sprint_matrices[l]
                to_invert = int_products[l]*(F0_p_mu^(p - d*l -1))
                acc = acc.apply_map(lambda x: x*to_invert^(-1))
                coefs = acc[0][d-g:]
                current_row += coefs

                for s in range(ell - 1, -1, -1): 
                    acc = acc * small_step_matrix[l][s].realize(p)
                    to_invert = small_step_denominators[l][s].realize(p)*(F0_p_mu^(d))
                    acc = acc.apply_map(lambda x: x*to_invert^(-1))
                    current_row += acc.list()
                
                # reversed: acc entries are ordered right-to-left in Harvey's
                # column convention, so the last entry corresponds to column 0
                A_f.append(list(reversed(current_row)))
            A_f_p = Matrix(A_f).apply_map(lambda x: Integer(x) % (p^(lam)))
            p_to_A_f[p] = A_f_p
    
    return p_to_A_f


def compute_A_f_avg_poly_from_curve(C, N, mu=2):
    F_coeffs_poly, _ = C.hyperelliptic_polynomials()
    d = C.degree()
    F_coeffs = [Integer(c) for c in F_coeffs_poly]
    return compute_A_F_l_avg_poly(F_coeffs, d, N, 1)


'''
AVG poly up to 50,000 took 59 seconds
Sqrt up to 50,000 took 1943 seconds (30 minutes) 
'''
        
    
import os as _os
if 'P8-avg-poly' in _os.path.basename(sys.argv[0]):
    start = timer()
    N = 199
    R = PolynomialRing(Integers(), 'x')
    x = R.gen()
    #f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
    f =  -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
    C = HyperellipticCurve(f)
    answer = compute_A_f_avg_poly_from_curve(C, N)
    for key, value in answer.items():
        print(key)
        print(value)

    time = timer() - start
    print(time)
