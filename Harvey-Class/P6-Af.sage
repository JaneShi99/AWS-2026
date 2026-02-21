import timeit as timeit_module 
timer = timeit_module.default_timer

mu = 3

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

def fast_product(lin_list):
    product_tree = build_product_tree(lin_list)
    return product_tree.val


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

def evaluation_tree(evaluation_poly, values, p):
    R.<k> = PolynomialRing(Zmod(p^mu))
    k = R.gen()
    leaves = [k-R(i) for i in values]
    prod_tree = build_product_tree(leaves)
    rem_tree, leaf_val_list = remainder_tree(evaluation_poly, prod_tree)
    return leaf_val_list

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

def compute_reusable_T_product(p, t, F_coeffs, m, d, R):
    # the fast matrix product T_(k+1)...T_(k+t) as a polynomial in k
    k = R.gen()
    T_matrices = [generate_T_matrix(p,k+i,d,m,F_coeffs) for i in range(1,t+1)]
    T_product = fast_product(T_matrices)

    return T_product

#compute T_(k0+1)...T_(k0+s)
def compute_T_product(p, k0, s, t, t_prime,F_coeffs, m, d, R, reusable_T):

    k = R.gen()

    # reusable_T is the fast matrix product T_(k+1)...T_(k+t) as a polynomial in k

    #reusable_T = compute_reusable_T_product(p, t, F_coeffs, m, d, R)

    values = [k0 + i*t for i in range(0, t)]
    nrows, ncols = reusable_T.nrows(), reusable_T.ncols()
    entry_evals = [evaluation_tree(reusable_T[i,j], values, p) for i in range(nrows) for j in range(ncols)]

    mat_vals = [
        matrix(R.base_ring(), nrows, ncols, [entry[k] for entry in entry_evals]) 
        for k in range(len(values))
    ]

    T_matrices_left = [generate_T_matrix(p, idx, d, m, F_coeffs) for idx in range(k0+t^2+1, k0+s+1)]

    total_mat_prod_t_squared = fast_product(mat_vals)
    if T_matrices_left:
        total_mat_prod_leftover = fast_product(T_matrices_left)
        total_mat_prod = total_mat_prod_t_squared * total_mat_prod_leftover
    else:
        total_mat_prod = total_mat_prod_t_squared


    # naive_prod = product([generate_T_matrix(p,i,d,m,F_coeffs) for i in range(k0+1,k0+s+1)])


    # print("total mat prod")
    # print(total_mat_prod)
    # print("naive prod")
    # print(product([generate_T_matrix(p,i,d,m,F_coeffs) for i in range(k0+1,k0+s+1)]))

    # assert total_mat_prod == naive_prod
    return total_mat_prod
 
# product of (k0+1, k0+2, ..., k0+s)
def compute_reusable_int(p, s, t, t_prime, R):
    k = R.gen()
    val_polys = [R(k + j) for j in range(1,t+1)]
    val_prod = fast_product(val_polys)
    return val_prod



def compute_int_products(p, k0, s, t, t_prime, R, reusable_int):
    t = floor(sqrt(s))
    t_prime = s - t^2
    k = R.gen()
    
    values = [R(k0 + j*t) for j in range(0, t)]
    evals = evaluation_tree(reusable_int, values, p)
    

    values_left = [R(j) for j in range(k0+t^2+1, k0+s+1)]
    final_product = product(evals) * product(values_left)
    
    return final_product


def divide_custom(x, y, p, R):
    divide_p_power = valuation(Integer(y), p)
    x_div_p = R(Integer(x)/(p^divide_p_power))
    y_div_p = R(Integer(y)/(p^divide_p_power))
    ans = x_div_p*(y_div_p^(-1))
    return ans


'''
Goal #1 is to just compute as is without reducing the extra O(g) factor
'''
    
def compute_A_f(F_coeffs, p):
    
    d = len(F_coeffs)-1
    g = (d-1) // 2
    m = (p-1)//2
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

    reusable_T = compute_reusable_T_product(p, t, F_coeffs, m, d, R)
    reusable_int = compute_reusable_int(p, s, t, t_prime, R)
    
    for i in range(0, g):
        # the following step is to compute the single step 
        # U_(ip-1)*T_(ip) = U_(ip)
        if i == 0:
            acc[0, -1] = R(F0^m)
            # print(acc)
        else:
            T_single = generate_T_matrix(p, p*i, d, m, F_coeffs)
            acc = acc * T_single
            to_invert = i*p*F0
            acc = acc.apply_map(lambda x: divide_custom(x, to_invert, p, Zmod(p^mu)))
            # print(acc)

        acc = acc * compute_T_product(p, i*p, s, t, t_prime, F_coeffs, m, d, R, reusable_T)

        # p-adic inversion is kind of broken here!
        
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
    

p = 59629079
R.<x> = PolynomialRing(Zmod(p^mu))
f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
discriminant(f)

start = timer()
compute_A_f(f.list(), p)
end = timer()
time = end - start
print(time)



#TODO: we should use the O(g) optimization