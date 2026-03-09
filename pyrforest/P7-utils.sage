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

# This can be generalized for R[P]/P^mu with mu being a variable.
class ZP2:
    """Element of R[P]/(P^2), stored as the 2x2 block matrix [[x, y], [0, x]].
    Represents x + y*P. Arithmetic is handled by native matrix operations."""
    __slots__ = ["mat"]

    def __init__(self, x, y=0):
        self.mat = matrix(ZZ, 2, 2, [x, y, 0, x])

    @property
    def x(self):
        return self.mat[0, 0]

    @property
    def y(self):
        return self.mat[0, 1]

    def __add__(self, other):
        result = ZP2.__new__(ZP2)
        result.mat = self.mat + other.mat
        return result

    def __sub__(self, other):
        result = ZP2.__new__(ZP2)
        result.mat = self.mat - other.mat
        return result

    def __mul__(self, other):
        result = ZP2.__new__(ZP2)
        result.mat = self.mat * other.mat
        return result

    def __neg__(self):
        result = ZP2.__new__(ZP2)
        result.mat = -self.mat
        return result

    def __mod__(self, n):
        result = ZP2.__new__(ZP2)
        result.mat = self.mat.apply_map(lambda v: v % n)
        return result

    def __repr__(self):
        return f"({self.x}) + ({self.y})*P"

    def __eq__(self, other):
        return self.mat == other.mat

    def scalar_mul(self, c):
        """Multiply by a scalar (element of the base ring)."""
        result = ZP2.__new__(ZP2)
        result.mat = c * self.mat
        return result

    def realize(self, p):
        """Evaluate at P=p to get a concrete value."""
        return self.x + self.y * p


# This can be generalized for R[P]/P^mu with mu being a variable.
class ZP2Matrix:
    """Matrix over R[P]/(P^2), stored as a single 2n x 2n block matrix
        [[m0, m1],
         [ 0, m0]]
    Represents M0 + M1*P. Arithmetic is handled by native matrix operations
    on the block matrix."""
    __slots__ = ["_mat", "_n_rows", "_n_cols"]

    def __init__(self, m0, m1):
        """Construct from the two component matrices m0, m1."""
        nr, nc = m0.nrows(), m0.ncols()
        self._n_rows = nr
        self._n_cols = nc
        self._mat = block_matrix([[m0, m1], [matrix(m0.base_ring(), nr, nc), m0]], subdivide=False)

    @staticmethod
    def _from_block(mat, nr, nc):
        """Internal: wrap an already-built 2n x 2n block matrix."""
        result = ZP2Matrix.__new__(ZP2Matrix)
        result._mat = mat
        result._n_rows = nr
        result._n_cols = nc
        return result

    @property
    def m0(self):
        """The constant part (top-left block)."""
        nr, nc = self._n_rows, self._n_cols
        return self._mat[:nr, :nc]

    @property
    def m1(self):
        """The coefficient of P (top-right block)."""
        nr, nc = self._n_rows, self._n_cols
        return self._mat[:nr, nc:]

    def __add__(self, other):
        return ZP2Matrix._from_block(self._mat + other._mat, self._n_rows, self._n_cols)

    def __sub__(self, other):
        return ZP2Matrix._from_block(self._mat - other._mat, self._n_rows, self._n_cols)

    def __mul__(self, other):
        # Block multiply: [[a0,a1],[0,a0]] * [[b0,b1],[0,b0]]
        # = [[a0*b0, a0*b1+a1*b0], [0, a0*b0]]  -- exactly the Z/P^2 rule
        return ZP2Matrix._from_block(self._mat * other._mat, self._n_rows, other._n_cols)

    def __neg__(self):
        return ZP2Matrix._from_block(-self._mat, self._n_rows, self._n_cols)

    def __mod__(self, n):
        return ZP2Matrix._from_block(self._mat.apply_map(lambda v: v % n), self._n_rows, self._n_cols)

    def scalar_mul(self, zp2_scalar):
        """Multiply every entry by a ZP2 scalar (x + y*P).
        Result: s0*[[m0,m1],[0,m0]] + s1*[[0,m0],[0,0]]
              = [[s0*m0, s0*m1+s1*m0], [0, s0*m0]]."""
        nr, nc = self._n_rows, self._n_cols
        base = self._mat.base_ring()
        Z = matrix(base, nr, nc)
        correction = block_matrix([[Z, self.m0], [Z, Z]], subdivide=False)
        return ZP2Matrix._from_block(
            zp2_scalar.x * self._mat + zp2_scalar.y * correction,
            nr, nc)

    def realize(self, p):
        """Evaluate at P=p to get a concrete Sage matrix."""
        return self.m0 + p * self.m1

    def __repr__(self):
        return f"ZP2Matrix (block form):\n  M0 = {self.m0}\n  M1 = {self.m1}"

    @staticmethod
    def zero(base_ring, nrows, ncols):
        """Create a zero ZP2Matrix."""
        return ZP2Matrix._from_block(matrix(base_ring, 2*nrows, 2*ncols), nrows, ncols)

    @staticmethod
    def from_array(base_ring, arr):
        """Convert an n x m list-of-lists of ZP2 elements into a ZP2Matrix.

        Usage:
            arr = [[ZP2(1,2), ZP2(3,4)],
                   [ZP2(5,6), ZP2(7,8)]]
            M = ZP2Matrix.from_array(ZZ, arr)
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
        return self._n_rows

    def ncols(self):
        return self._n_cols