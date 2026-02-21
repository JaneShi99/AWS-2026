import timeit as timeit_module
class TreeNode:
    def __init__(self, val=0, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right
    

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

def compute_Q_k(p, t):
    R.<k> = PolynomialRing(Zmod(p^2))
    lin_list = [k+i for i in range(1, t+1)]
    tree = build_product_tree(lin_list)
    return tree, tree.val

tree, res = compute_Q_k(11, 16)

def remainder_tree(root_poly, prod_tree_root):
    leaf_val_list = []
    def recursive_helper(poly, prod_tree): 
        if prod_tree is None: return None
        root_val = poly % prod_tree.val 
        left_rem = recursive_helper(root_val, prod_tree.left)
        right_rem = recursive_helper(root_val, prod_tree.right)
        if left_rem is None and right_rem is None:
            leaf_val_list.append(root_val)
        return TreeNode(root_val, left_rem, right_rem)
    
    rem_tree = recursive_helper(root_poly, prod_tree_root)
    return rem_tree, leaf_val_list

def compute_factorial(p):
    s = p-1
    t = floor(sqrt(s))
    t_prime = s - t^2


    R.<k> = PolynomialRing(Zmod(p^2))

    tree, Q_k = compute_Q_k(p, t)
    prod_tree = build_product_tree([k-i*t for i in range(0, t)])

    rem_tree, leaf_val_list = remainder_tree(Q_k, prod_tree)
    return product(leaf_val_list) * product(i for i in range(t^2+1, p))

def compute_factorial_naively(p):
    accumulator = 1
    for i in range(1, p):
        accumulator = (accumulator * i) % p^2
    return accumulator

def test_two_methods():
    exponents = [i for i in range(10, 30)]
    timer = timeit_module.default_timer

    for exp in exponents:
        p = next_prime(2^exp)
        print(f"p = {p} (2^{exp})")

        start1 = timer()
        res1 = compute_factorial(p)
        end1 = timer()
        t1 = end1 - start1

        start2 = timer()
        res2 = compute_factorial_naively(p)
        end2 = timer()
        t2 = end2 - start2
        
        print(f"  Tree method result:   {res1}")
        print(f"  Direct method result: {res2}")
        print(f"  Match: {res1 == res2}")

        print(f"  Tree method:   {t1:.6f}s")

        print(f"  Direct method: {t2:.6f}s")
        print()

#test_two_methods()