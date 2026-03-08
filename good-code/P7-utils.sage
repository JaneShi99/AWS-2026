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

# This can be genealized for R[P]/P^mu with mu being a variable.
class ZP2:
    """Element of R[P]/(P^2), stored as (x, y) representing x + y*P."""
    __slots__ = ["x", "y"]

    def __init__(self, x, y=0):
        self.x = x
        self.y = y

    def __add__(self, other):
        return ZP2(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return ZP2(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return ZP2(self.x * other.x, self.x * other.y + self.y * other.x)

    def __neg__(self):
        return ZP2(-self.x, -self.y)

    def __mod__(self, n):
        return ZP2(self.x % n, self.y % n)

    def __repr__(self):
        return f"({self.x}) + ({self.y})*P"

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def scalar_mul(self, c):
        """Multiply by a scalar (element of the base ring)."""
        return ZP2(c * self.x, c * self.y)

    def realize(self, p):
        """Evaluate at P=p to get a concrete value."""
        return self.x + self.y * p


# This can be genealized for R[P]/P^mu with mu being a variable.
class ZP2Matrix:
    """Matrix over R[P]/(P^2), stored as two Sage matrices (M0, M1).
    Represents M0 + M1*P. Uses Sage's fast matrix arithmetic internally."""
    __slots__ = ["m0", "m1"]

    def __init__(self, m0, m1):
        self.m0 = m0  # constant part (Sage matrix over base ring)
        self.m1 = m1  # coefficient of P (Sage matrix over base ring)

    def __add__(self, other):
        return ZP2Matrix(self.m0 + other.m0, self.m1 + other.m1)

    def __sub__(self, other):
        return ZP2Matrix(self.m0 - other.m0, self.m1 - other.m1)

    def __mul__(self, other):
        # (A0 + A1*P)(B0 + B1*P) = A0*B0 + (A0*B1 + A1*B0)*P
        return ZP2Matrix(
            self.m0 * other.m0,
            self.m0 * other.m1 + self.m1 * other.m0
        )

    def __neg__(self):
        return ZP2Matrix(-self.m0, -self.m1)

    def __mod__(self, n):
        return ZP2Matrix(self.m0.apply_map(lambda x: x % n), self.m1.apply_map(lambda x: x % n))

    def scalar_mul(self, zp2_scalar):
        """Multiply every entry by a ZP2 scalar (x + y*P)."""
        return ZP2Matrix(
            zp2_scalar.x * self.m0,
            zp2_scalar.x * self.m1 + zp2_scalar.y * self.m0
        )

    def realize(self, p):
        """Evaluate at P=p to get a concrete Sage matrix."""
        return self.m0 + p * self.m1

    def __repr__(self):
        return f"ZP2Matrix:\n  M0 = {self.m0}\n  M1 = {self.m1}"

    @staticmethod
    def zero(base_ring, nrows, ncols):
        """Create a zero ZP2Matrix."""
        z = matrix(base_ring, nrows, ncols)
        return ZP2Matrix(z, copy(z))

    @staticmethod
    def from_array(base_ring, arr):
        """Convert an n×n list-of-lists of ZP2 elements into a ZP2Matrix.
        
        Usage:
            arr = [[ZP2(1,2), ZP2(3,4)],
                   [ZP2(5,6), ZP2(7,8)]]
            M = ZP2Matrix.from_array(Zmod(p^3), arr)
        """
        nrows = len(arr)
        ncols = len(arr[0])
        m0 = matrix(base_ring, nrows, ncols, [arr[i][j].x for i in range(nrows) for j in range(ncols)])
        m1 = matrix(base_ring, nrows, ncols, [arr[i][j].y for i in range(nrows) for j in range(ncols)])
        return ZP2Matrix(m0, m1)

    @staticmethod
    def identity(base_ring, n):
        """Create an identity ZP2Matrix."""
        return ZP2Matrix(matrix.identity(base_ring, n), matrix(base_ring, n, n))

    def nrows(self):
        return self.m0.nrows()

    def ncols(self):
        return self.m0.ncols()