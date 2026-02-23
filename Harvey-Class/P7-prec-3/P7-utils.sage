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


class ZP3:
    """Element of R[P]/(P^3), stored as (x, y, z) representing x + y*P + z*P^2."""
    __slots__ = ["x", "y", "z"]

    def __init__(self, x, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return ZP3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return ZP3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        # (a + b*P + c*P^2)(d + e*P + f*P^2) = ad + (ae+bd)*P + (af+be+cd)*P^2  mod P^3
        return ZP3(
            self.x * other.x,
            self.x * other.y + self.y * other.x,
            self.x * other.z + self.y * other.y + self.z * other.x
        )

    def __neg__(self):
        return ZP3(-self.x, -self.y, -self.z)

    def __mod__(self, n):
        return ZP3(self.x % n, self.y % n, self.z % n)

    def __repr__(self):
        return f"({self.x}) + ({self.y})*P + ({self.z})*P^2"

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def scalar_mul(self, c):
        """Multiply by a scalar (element of the base ring)."""
        return ZP3(c * self.x, c * self.y, c * self.z)

    def realize(self, p):
        """Evaluate at P=p to get a concrete value."""
        return self.x + self.y * p + self.z * p^2


class ZP3Matrix:
    """Matrix over R[P]/(P^3), stored as three Sage matrices (M0, M1, M2).
    Represents M0 + M1*P + M2*P^2. Uses Sage's fast matrix arithmetic internally."""
    __slots__ = ["m0", "m1", "m2"]

    def __init__(self, m0, m1, m2):
        self.m0 = m0  # constant part (Sage matrix over base ring)
        self.m1 = m1  # coefficient of P (Sage matrix over base ring)
        self.m2 = m2  # coefficient of P^2 (Sage matrix over base ring)

    def __add__(self, other):
        return ZP3Matrix(self.m0 + other.m0, self.m1 + other.m1, self.m2 + other.m2)

    def __sub__(self, other):
        return ZP3Matrix(self.m0 - other.m0, self.m1 - other.m1, self.m2 - other.m2)

    def __mul__(self, other):
        # (A0 + A1*P + A2*P^2)(B0 + B1*P + B2*P^2)
        # = A0*B0 + (A0*B1 + A1*B0)*P + (A0*B2 + A1*B1 + A2*B0)*P^2  mod P^3
        return ZP3Matrix(
            self.m0 * other.m0,
            self.m0 * other.m1 + self.m1 * other.m0,
            self.m0 * other.m2 + self.m1 * other.m1 + self.m2 * other.m0
        )

    def __neg__(self):
        return ZP3Matrix(-self.m0, -self.m1, -self.m2)

    def __mod__(self, n):
        return ZP3Matrix(
            self.m0.apply_map(lambda x: x % n),
            self.m1.apply_map(lambda x: x % n),
            self.m2.apply_map(lambda x: x % n)
        )

    def scalar_mul(self, zp3_scalar):
        """Multiply every entry by a ZP3 scalar (x + y*P + z*P^2)."""
        return ZP3Matrix(
            zp3_scalar.x * self.m0,
            zp3_scalar.x * self.m1 + zp3_scalar.y * self.m0,
            zp3_scalar.x * self.m2 + zp3_scalar.y * self.m1 + zp3_scalar.z * self.m0
        )

    def realize(self, p):
        """Evaluate at P=p to get a concrete Sage matrix."""
        return self.m0 + p * self.m1 + p^2 * self.m2

    def __repr__(self):
        return f"ZP3Matrix:\n  M0 = {self.m0}\n  M1 = {self.m1}\n  M2 = {self.m2}"

    @staticmethod
    def zero(base_ring, nrows, ncols):
        """Create a zero ZP3Matrix."""
        z = matrix(base_ring, nrows, ncols)
        return ZP3Matrix(z, copy(z), copy(z))

    @staticmethod
    def from_array(base_ring, arr):
        """Convert an n×n list-of-lists of ZP3 elements into a ZP3Matrix.
        
        Usage:
            arr = [[ZP3(1,2,3), ZP3(4,5,6)],
                   [ZP3(7,8,9), ZP3(10,11,12)]]
            M = ZP3Matrix.from_array(ZZ, arr)
        """
        nrows = len(arr)
        ncols = len(arr[0])
        m0 = matrix(base_ring, nrows, ncols, [arr[i][j].x for i in range(nrows) for j in range(ncols)])
        m1 = matrix(base_ring, nrows, ncols, [arr[i][j].y for i in range(nrows) for j in range(ncols)])
        m2 = matrix(base_ring, nrows, ncols, [arr[i][j].z for i in range(nrows) for j in range(ncols)])
        return ZP3Matrix(m0, m1, m2)

    @staticmethod
    def identity(base_ring, n):
        """Create an identity ZP3Matrix."""
        return ZP3Matrix(matrix.identity(base_ring, n), matrix(base_ring, n, n), matrix(base_ring, n, n))

    def nrows(self):
        return self.m0.nrows()

    def ncols(self):
        return self.m0.ncols()