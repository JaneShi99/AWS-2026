import timeit as timeit_module 
timer = timeit_module.default_timer

mu = 3
# really seem like we want mu = 3. mu=2 doesn't seem to match for say 101.
# ?? TODO: maybe ask David about this.
class TreeNode:
    def __init__(self, val=0, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right



'''
Attempts to print a tree, but is terrible
'''

def print_tree(root):
    if root is None: 
        return 
    
    q = [root]
    while q:
        next_q = []
        for node in q:
            print(node.val, end=", ")
            if node.left: next_q.append(node.left)
            if node.right: next_q.append(node.right)
        print()
        q = next_q


'''
Given a list of objects, build the product tree. 
The product is not type dependent.
'''
def build_product_tree(lin_list):
    if not lin_list: return None
    
    def recursive_helper(a, b):
        if a == b: return None
        elif a == b-1: return TreeNode(lin_list[a])
        else:
            m = a + (b-a)//2 
            left = recursive_helper(a, m)
            right = recursive_helper(m, b)
            left_val = left.val if left else 1
            right_val = right.val if right else 1
            return TreeNode(left_val * right_val, left, right)
    
    return recursive_helper(0, len(lin_list))


'''
Given a list of objects, compute the product of the objects.
'''
def fast_product(lin_list):
    product_tree = build_product_tree(lin_list)
    return product_tree.val

'''
Given a polynomial (as the root), and a product tree coming from a list of linear polynomials,
compute the modular evaluation top-down.
'''
def remainder_tree(root_element, prod_tree_root):
    leaf_val_list = []
    def recursive_helper(element, prod_tree): 
        if prod_tree is None: return None
        root_val = element % prod_tree.val 
        left_rem = recursive_helper(root_val, prod_tree.left)
        right_rem = recursive_helper(root_val, prod_tree.right)
        if left_rem is None and right_rem is None:
            leaf_val_list.append(root_val)
        return TreeNode(root_val, left_rem, right_rem)
    
    rem_tree = recursive_helper(root_element, prod_tree_root)
    return rem_tree, leaf_val_list

'''
Given a list of values, and a polynomial, using the remainder tree, 
compute the polynomial evaluated at each of the values.
'''
def evaluation_tree(evaluation_poly, values, p):
    R.<k> = PolynomialRing(Zmod(p^mu))
    k = R.gen()
    leaves = [k-R(i) for i in values]
    prod_tree = build_product_tree(leaves)
    rem_tree, leaf_val_list = remainder_tree(evaluation_poly, prod_tree)
    return leaf_val_list

'''
Compute the single matrix T_k
'''
def generate_T_matrix(p, k_val, d, m, F_coeffs):
    R.<k> = PolynomialRing(Zmod(p^mu))
    sub = R(k_val) * F_coeffs[0]
    M = matrix(R, d, d)
    for i in range(1,d):
        M[i,i-1] = sub
    
    for i in range(d):
        j = d - i
        M[i, d-1] = (R(j*(m+1)) - R(k_val)) * F_coeffs[j]

    return M

'''
Compute the product of matrices T_(k+1)...T_(k+t) as a matrix whose entries are polynomials in k.
Using fast matrix product.
'''
def compute_reusable_T_product(p, t, F_coeffs, m, d, R):
    k = R.gen()
    T_matrices = [generate_T_matrix(p,k+i,d,m,F_coeffs) for i in range(1,t+1)]
    T_product = fast_product(T_matrices)

    return T_product

'''
reusable_T is the output of the above function, which is only evaluated once in the algorithm.
In this function, we compute the product of matrices T_(k0+1)...T_(k0+s)
by plugging in those values and using the evaluation tree.
'''
def compute_T_product(p, k0, s, t, t_prime,F_coeffs, m, d, R, reusable_T):

    k = R.gen()

    values = [k0 + i*t for i in range(0, t)]
    nrows, ncols = reusable_T.nrows(), reusable_T.ncols()

    #using evaluation tree to compute the matrix product T_(k0+1)...T_(k0+t^2)
    entry_evals = [evaluation_tree(reusable_T[i,j], values, p) for i in range(nrows) for j in range(ncols)]
    mat_vals = [
        matrix(R.base_ring(), nrows, ncols, [entry[k] for entry in entry_evals]) 
        for k in range(len(values))
    ]
    total_mat_prod_t_squared = fast_product(mat_vals)

    #compute the product of the "leftover" matrices T_(k0+t^2+1)...T_(k0+s)
    T_matrices_left = [generate_T_matrix(p, idx, d, m, F_coeffs) for idx in range(k0+t^2+1, k0+s+1)]

    if T_matrices_left:
        total_mat_prod_leftover = fast_product(T_matrices_left)
        total_mat_prod = total_mat_prod_t_squared * total_mat_prod_leftover
    else:
        total_mat_prod = total_mat_prod_t_squared

    return total_mat_prod
 
# product of (k+1, k+2, ... k+t) as a polynomial in k
def compute_reusable_int(p, s, t, t_prime, R):
    k = R.gen()
    val_polys = [R(k + j) for j in range(1,t+1)]
    val_prod = fast_product(val_polys)
    return val_prod

# product of (k0+1, k0+2, ... k0+s) as a polynomial in k
def compute_int_products(p, k0, s, t, t_prime, R, reusable_int):
    t = floor(sqrt(s))
    t_prime = s - t^2
    k = R.gen()
    
    #compute the product of (k0+1, k0+2, ... k0+t^2) using reusable_int
    values = [R(k0 + j*t) for j in range(0, t)]
    evals = evaluation_tree(reusable_int, values, p)
    
    #compute the product of (k0+t^2+1, k0+t^2+2, ... k0+s)
    values_left = [R(j) for j in range(k0+t^2+1, k0+s+1)]
    final_product = product(evals) * product(values_left)
    
    return final_product

# custom division in Z/p^mu Z
def divide_custom(x, y, p, R):
    divide_p_power = valuation(Integer(y), p)
    x_div_p = R(Integer(x)/(p^divide_p_power))
    y_div_p = R(Integer(y)/(p^divide_p_power))
    ans = R(x_div_p*(y_div_p^(-1)))
    return ans


'''
Goal #1 is to just compute as is without reducing the extra O(g) factor
'''
    
def compute_A_f(F_coeffs, p):
    
    d = len(F_coeffs)-1
    g = (d-1) // 2
    m = (p-1) // 2
    R.<k> = PolynomialRing(Zmod(p^mu))
    k = R.gen()
    F_coeffs = [R(c) for c in F_coeffs]
    F0 = R(F_coeffs[0])

    acc = Matrix(R, 1, d)
    A_f = []

    F0_p = F_coeffs[0]^(p-1)

    s = p-1
    t = floor(sqrt(s))
    t_prime = s - t^2

    # the following are reusable objects for the sprint matrix products and 
    # for the constant products
    reusable_T = compute_reusable_T_product(p, t, F_coeffs, m, d, R)
    reusable_int = compute_reusable_int(p, s, t, t_prime, R)
    
    for i in range(0, g):
        # the following step is to compute the single step 
        # U_(ip) = U_(ip-1)*T_(ip)*(constants)

        if i == 0:
            acc[0, -1] = R(F0^m)
            # print(acc)
        else:
            T_single = generate_T_matrix(p, p*i, d, m, F_coeffs)
            acc = acc * T_single
            to_invert = i*p*F0
            acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, Zmod(p^mu)))
            # print(acc)

        # the following step is to compute the sprint step 
        # U_((i+1)*p) = U_(ip)*T_(ip+1)*T_(ip+2)*...*T_(ip+p-1) * constants
        acc = acc * compute_T_product(p, i*p, s, t, t_prime, F_coeffs, m, d, R, reusable_T)
        to_invert = (Zmod(p^mu)(compute_int_products(p, i*p, s, t, t_prime, R, reusable_int)*F0_p))
        acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, Zmod(p^mu)))
        
        A_f.append(list(acc[0][(-g):]))


    print("computed A_f")
    A_mod = Matrix(A_f).apply_map(lambda x: Integer(x) % p)
    print(A_mod)
    return A_f


'''
Goal #2 is to just compute with the extra O(g) factor
'''
def compute_A_f_fast(F_coeffs, p):
    
    d = len(F_coeffs)-1
    g = (d-1) // 2
    m = (p-1) // 2
    R.<k> = PolynomialRing(Zmod(p^mu))
    k = R.gen()
    F_coeffs = [R(c) for c in F_coeffs]
    F0 = R(F_coeffs[0])

    acc = Matrix(R, 1, d)
    A_f = []

    F0_p = F_coeffs[0]^(p-1)

    s = p-1
    t = floor(sqrt(s))
    t_prime = s - t^2

    # the following are reusable objects for the sprint matrix products and 
    # for the constant products
    reusable_T = compute_reusable_T_product(p, t, F_coeffs, m, d, R)
    reusable_int = compute_reusable_int(p, s, t, t_prime, R)

    R_0 = compute_T_product(p, 0, s, t, t_prime, F_coeffs, m, d, R, reusable_T)
    R_p = compute_T_product(p, p, s, t, t_prime, F_coeffs, m, d, R, reusable_T)

    R0_minus_Rp = R_p - R_0
    R_1 = R0_minus_Rp.apply_map(lambda x: R(Integer(x)/p))
    sprint_matrices = [R_0 + (i*p)*R_1 for i in range(0, g)]

    print(p)
    print("sprint")
    for mat in sprint_matrices:
        print(mat.apply_map(lambda x: Integer(x) % (p^2)))
        print("\n")

    int_0 = compute_int_products(p, 0, s, t, t_prime, R, reusable_int)
    int_p = compute_int_products(p, p, s, t, t_prime, R, reusable_int)
    int_0_minus_int_p = int_p - int_0
    int_1 = R(Integer(int_0_minus_int_p)/p)

    int_products = [int_0 + (i*p)*int_1 for i in range(0, g)]

    #print("int")
    #print(int_products)


    for i in range(0, g):
        # the following step is to compute the single step 
        # U_(ip) = U_(ip-1)*T_(ip)*(constants)

        if i == 0:
            acc[0, -1] = R(F0^m)
            # print(acc)
        else:
            T_single = generate_T_matrix(p, p*i, d, m, F_coeffs)
            acc = acc * T_single
            to_invert = i*p*F0
            #print(to_invert)
            acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, Zmod(p^mu)))
            # print(acc)

        # the following step is to compute the sprint step 
        # U_((i+1)*p) = U_(ip)*T_(ip+1)*T_(ip+2)*...*T_(ip+p-1) * constants
        acc = acc * sprint_matrices[i]
        to_invert = int_products[i]*F0^(p-1)
        acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, Zmod(p^mu)))
        
        A_f.append(list(acc[0][(-g):]))


    A_mod = Matrix(A_f).apply_map(lambda x: Integer(x) % (p))
    
    return A_mod
    
def exp1():
    p = 4999
    R.<x> = PolynomialRing(Zmod(p^mu))
    f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
    discriminant(f)


    start = timer()
    ans = compute_A_f(f.list(), p)
    print(ans)
    end = timer()
    time = end - start
    print(time)

    mu = 2

    start = timer()
    ans2 = compute_A_f_fast(f.list(), p)
    print(ans2)
    end = timer()
    time = end - start
    print(time)

'''
[ 8  2  8 91  5]
[66 94 77 31  6]
[36 65 56  9 16]
[16  8 87 67 75]
[19 56 56 28 78]
'''
def exp2():
    N = 102

    start = timer()
    for p in range(5, N):
        if is_prime(p):
            print(p)
            mu = 2
            R.<x> = PolynomialRing(Zmod(p^mu))
            #f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
            f =  -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
            discriminant(f)

            ans = compute_A_f_fast(f.list(), p)
            print(ans)
    
    end = timer()
    time = end - start
    print(time)

'''
AVG poly up to 50,000 took 59 seconds
Sqrt up to 50,000 took 1943 seconds (30 minutes) 
'''

# exp2()



def exp3():
    N = 6299
    start = timer()
    for p in [6299]:
        if is_prime(p):
            print(p)
            R.<x> = PolynomialRing(Zmod(p^mu))
            #f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
            f =  -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
            discriminant(f)

            ans = compute_A_f_fast(f.list(), p)
            print(ans)
    
    end = timer()
    time = end - start
    print(time)

'''
AVG poly up to 50,000 took 59 seconds
Sqrt up to 50,000 took 1943 seconds (30 minutes) 
'''

exp3()