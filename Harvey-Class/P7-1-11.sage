load("P6-tree.sage")

import timeit as timeit_module

timer = timeit_module.default_timer

'''
The following implements ``accumulating remainder tree``
'''
def remainder_tree_builder(value_tree, modulus_tree):
    leave_val_list = []
    def recursive_helper(value_subtree, modulus_subtree, higher_value, left_sub, left_nbd_val):
        if value_subtree is None: return None

        root_val = (higher_value if left_sub else higher_value * left_nbd_val) % modulus_subtree.val 

        left_rem = recursive_helper(value_subtree.left, modulus_subtree.left, root_val, True, -1)

        right_rem = recursive_helper(value_subtree.right, modulus_subtree.right, root_val, False, value_subtree.left.val if value_subtree.left else -1)

        if left_rem is None and right_rem is None:
            leave_val_list.append(root_val)

        return TreeNode(root_val, left_rem, right_rem)

    return recursive_helper(value_tree, modulus_tree, 1, True, -1), leave_val_list


def wilson(p):
    value_tree = build_product_tree([k for k in range(1,p+1)])
    modulus_tree = build_product_tree([k^2 if is_prime(k) else 1 for k in range(1,p+1)])
    
    rem_tree, leaf_val_list = remainder_tree_builder(value_tree, modulus_tree)
    
    return leaf_val_list

def test_two_methods():
    n = 50000
    start1 = timer()
    ans = (wilson(n))
    end1 = timer()
    print(f"fast method: {end1-start1}")


    start2 = timer()
    for i in range(1, n+1):
        if is_prime(i):
            compute_factorial(i)
    end2 = timer()
    print(f"slow method: {end2-start2}")
test_two_methods()


'''
n=50000
fast method: 0.6834103129804134
slow method: 10.797543739434332
'''