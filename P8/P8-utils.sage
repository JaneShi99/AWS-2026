import timeit as timeit_module
timer = timeit_module.default_timer

class TreeNode:
    def __init__(self, val=0, left=None, right=None):
        self.val = val
        self.left = left
        self.right = right


def build_product_tree(lin_list):
    if not lin_list: return None
    
    def recursive_helper(a, b):
        if a == b: return None
        elif a == b-_sage_const_1 : return TreeNode(lin_list[a])
        else:
            m = a + (b-a)//_sage_const_2  
            left = recursive_helper(a, m)
            right = recursive_helper(m, b)
            left_val = left.val if left else _sage_const_1 
            right_val = right.val if right else _sage_const_1 
            return TreeNode(left_val * right_val, left, right)
    
    return recursive_helper(_sage_const_0 , len(lin_list))

'''
The following implements ``accumulating remainder tree``
'''
def remainder_tree_builder(value_tree, modulus_tree, identity = 1):
    leave_val_list = []
    def recursive_helper(value_subtree, modulus_subtree, higher_value, left_sub, left_nbd_val):
        if value_subtree is None: return None

        root_val = (higher_value if left_sub else higher_value * left_nbd_val) % modulus_subtree.val 

        left_rem = recursive_helper(value_subtree.left, modulus_subtree.left, root_val, True, -1)

        right_rem = recursive_helper(value_subtree.right, modulus_subtree.right, root_val, False, value_subtree.left.val if value_subtree.left else -1)

        if left_rem is None and right_rem is None:
            leave_val_list.append(root_val)

        return TreeNode(root_val, left_rem, right_rem)

    return recursive_helper(value_tree, modulus_tree, identity, True, -1), leave_val_list

class ZPmu:
    """Element of R[P]/(P^mu), stored as a list of mu coefficients.
    Represents coeffs[0] + coeffs[1]*P + ... + coeffs[mu-1]*P^(mu-1)."""
    __slots__ = ["coeffs", "mu"]

    def __init__(self, coeffs, mu):
        """coeffs: list of length <= mu. Padded with 0 to length mu."""
        self.mu = mu
        self.coeffs = list(coeffs) + [0] * (mu - len(coeffs))

    def __add__(self, other):
        assert self.mu == other.mu
        return ZPmu([self.coeffs[i] + other.coeffs[i] for i in range(self.mu)], self.mu)

    def __sub__(self, other):
        assert self.mu == other.mu
        return ZPmu([self.coeffs[i] - other.coeffs[i] for i in range(self.mu)], self.mu)

    def __mul__(self, other):
        """Truncated convolution mod P^mu."""
        assert self.mu == other.mu
        mu = self.mu
        result = [0] * mu
        for k in range(mu):
            s = 0
            for i in range(k + 1):
                s += self.coeffs[i] * other.coeffs[k - i]
            result[k] = s
        return ZPmu(result, mu)

    def __neg__(self):
        return ZPmu([-c for c in self.coeffs], self.mu)

    def __mod__(self, n):
        return ZPmu([c % n for c in self.coeffs], self.mu)

    def __repr__(self):
        terms = [f"({self.coeffs[i]})*P^{i}" for i in range(self.mu)]
        return " + ".join(terms)

    def __eq__(self, other):
        return self.mu == other.mu and self.coeffs == other.coeffs

    def scalar_mul(self, c):
        """Multiply by a scalar (element of the base ring)."""
        return ZPmu([c * ci for ci in self.coeffs], self.mu)

    def realize(self, p):
        """Evaluate at P=p to get a concrete value."""
        return sum(self.coeffs[i] * p^i for i in range(self.mu))


class ZPmuMatrix:
    """Matrix over R[P]/(P^mu), stored as a list of mu Sage matrices.
    Represents mats[0] + mats[1]*P + ... + mats[mu-1]*P^(mu-1).
    Uses Sage's fast matrix arithmetic internally."""
    __slots__ = ["mats", "mu"]

    def __init__(self, mats, mu=None):
        """mats: list of Sage matrices of length <= mu. Padded with zero matrices to length mu.
        If mu is None, infer from len(mats)."""
        if mu is None:
            mu = len(mats)
        self.mu = mu
        if len(mats) < mu:
            z = matrix(mats[0].base_ring(), mats[0].nrows(), mats[0].ncols())
            self.mats = list(mats) + [copy(z) for _ in range(mu - len(mats))]
        else:
            self.mats = list(mats)

    def __add__(self, other):
        assert self.mu == other.mu
        return ZPmuMatrix([self.mats[i] + other.mats[i] for i in range(self.mu)], self.mu)

    def __sub__(self, other):
        assert self.mu == other.mu
        return ZPmuMatrix([self.mats[i] - other.mats[i] for i in range(self.mu)], self.mu)

    def __mul__(self, other):
        """Truncated convolution of matrices mod P^mu."""
        assert self.mu == other.mu
        mu = self.mu
        result = []
        for k in range(mu):
            s = self.mats[0] * other.mats[k] if k == 0 else sum(
                self.mats[i] * other.mats[k - i] for i in range(k + 1)
            )
            result.append(s)
        return ZPmuMatrix(result, mu)

    def __neg__(self):
        return ZPmuMatrix([-m for m in self.mats], self.mu)

    def __mod__(self, n):
        return ZPmuMatrix([m.apply_map(lambda x: x % n) for m in self.mats], self.mu)

    def scalar_mul(self, zpmu_scalar):
        """Multiply every entry by a ZPmu scalar via truncated convolution."""
        assert self.mu == zpmu_scalar.mu
        mu = self.mu
        result = []
        for k in range(mu):
            s = zpmu_scalar.coeffs[0] * self.mats[k] if k == 0 else sum(
                zpmu_scalar.coeffs[i] * self.mats[k - i] for i in range(k + 1)
            )
            result.append(s)
        return ZPmuMatrix(result, mu)

    def realize(self, p):
        """Evaluate at P=p to get a concrete Sage matrix."""
        return sum(self.mats[i] * p^i for i in range(self.mu))

    def __repr__(self):
        parts = [f"  M{i} = {self.mats[i]}" for i in range(self.mu)]
        return "ZPmuMatrix:\n" + "\n".join(parts)

    @staticmethod
    def zero(base_ring, nrows, ncols, mu):
        """Create a zero ZPmuMatrix."""
        return ZPmuMatrix([matrix(base_ring, nrows, ncols) for _ in range(mu)], mu)

    @staticmethod
    def from_array(base_ring, arr, mu):
        """Convert an n×n list-of-lists of ZPmu elements into a ZPmuMatrix.
        
        Usage:
            arr = [[ZPmu([1,2], 2), ZPmu([3,4], 2)],
                   [ZPmu([5,6], 2), ZPmu([7,8], 2)]]
            M = ZPmuMatrix.from_array(ZZ, arr, 2)
        """
        nrows = len(arr)
        ncols = len(arr[0])
        inferred_mu = arr[0][0].mu
        assert inferred_mu == mu
        mats = []
        for idx in range(mu):
            m = matrix(base_ring, nrows, ncols,
                       [arr[i][j].coeffs[idx] for i in range(nrows) for j in range(ncols)])
            mats.append(m)
        return ZPmuMatrix(mats, mu)

    @staticmethod
    def identity(base_ring, n, mu):
        """Create an identity ZPmuMatrix."""
        mats = [matrix.identity(base_ring, n)] + [matrix(base_ring, n, n) for _ in range(mu - 1)]
        return ZPmuMatrix(mats, mu)

    def nrows(self):
        return self.mats[0].nrows()

    def ncols(self):
        return self.mats[0].ncols()