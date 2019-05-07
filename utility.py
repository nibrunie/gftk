# -*- coding: utf-8 -*-

import random


def clm(a, b):
    """ generic carry-less multiply """
    c = b
    acc = 0
    i = 0
    while c >> i!= 0:
        if (c >> i) & 1:
            acc ^= a << i
        i += 1
    return acc

def gf_mod(a, r, m=64, n=32):
    """ generic galois-field modulo return a [X^n+r], deg(a) < m """
    assert m >= n
    mod = (r << (m - n)) ^ (1 << m)
    mask = (2^n - 1)
    rem = a
    for i in range(m,n-1,-1):
        if rem & (1 << i):
            rem ^= mod
        mod = mod >> 1
    return rem

def gf_mul(a, b, r,m0=32,m1=32,n=32):
    """ Galois field multiplication
        @p a is at most m0-bit wide
        @p b is at most m1-bit wide
        @p r is at most n-bit wide
        return clm(a, b) [X^n+r]"""
    return gf_mod(clm(a, b), r, (m0+m1), n)


def value_to_bit_matrix(a, s=(8, 8)):
    """ convert a from n-bit wide value to
        s[0] x s[1] bit matrix"""
    bit_matrix = [[0] * s[1] for i in range(s[0])]
    for i in range(s[0]):
        for j in range(s[1]):
            index = j + i * s[0]
            bit = (a & (1 << index)) >> index
            bit_matrix[i][j] = bit
    return bit_matrix

def value_to_bit_vector(v, n=8):
    """ Convert value @p v to bit Vector of size @p n """
    return Vector([(v >> i) & 1 for i in range(n)])


def gf_mul_by_cst(v, cst, r, m0=8, m1=8, n=8):
    """ generate BitMatrix-based multiplication of @p v . @p cst """
    # building bit matrix
    bm = BitMatrix.from_value_list([gf_mul(cst, 1 << i, r, m0, n, n) for i in range(n)])
    print("bm=\n{}".format(str(bm)))
    value_vector = value_to_bit_vector(v, m0)
    result_vector = value_vector * bm
    return result_vector.value



def bit_matrix_multiply(a, b):
    """ Matrix multiplication between two matrices """
    assert a.col == b.row
    result = BitMatrix(0, a.row, b.col)
    for i in range(a.row):
        for j in range(b.col):
            acc = 0
            for k in range(a.col):
                acc ^= a[i][k] * b[k][j]
            result[i][j] = acc
    return result

class BitMatrix:
    """ Bit Matrix structure """
    def __init__(self, value, row=8, col=8):
        self.matrix = value_to_bit_matrix(value, (row, col))
        self.row = row
        self.col = col

    def __getitem__(self, row_id):
        """ select a matrix row """
        return self.matrix[row_id]

    def __setitem__(self, row_id, row_value):
        """ set a matrix row (by copying value from @p row_value) """
        self.matrix[row_id] = [v for v in row_value]

    def __mul__(self, b):
        """ implicit matrix multiply """
        return bit_matrix_multiply(self, b)

    @staticmethod
    def from_value_list(value_list, elt_size=8):
        """ transform a list of vector into a matrix,
            each vector is stored as a row """
        row = len(value_list)
        result = BitMatrix(0, row, elt_size)
        for i in range(row):
            result[i] = value_to_bit_vector(value_list[i], elt_size)
        return result

    @staticmethod
    def identity(row_col):
        """ generate identity matrix of size @p row_col x row_col """
        result = BitMatrix(0, row_col, row_col)
        for i in range(row_col):
            result[i][i] = 1
        return result

    @property
    def transpose(self):
        """ generate the transpose of matrix @p self """
        result = BitMatrix(0, self.col, self.row)
        for i in range(self.row):
            for j in range(self.col):
                result[i][j] = self[j][i]
        return result

    @property
    def scalar_value(self):
        """ convert bit matrix to scalar value """
        acc = 0
        for i in range(self.row):
            for j in range(self.col):
                index = j + i * self.row
                acc ^= 1 << index
        return acc

    def vector_row(self, index):
        """ return @p index-th row of self matrix """
        return Vector(self.matrix[index])

    def __str__(self):
        """ Convert matrix to string representation """
        return "\n".join(" ".join("{}".format(self.matrix[i][j]) for j in range(self.col)) for i in range(self.row))

def build_gen_vector(k, elt, red_pol=0x1d, elt_size=8, red_pol_size=8):
    """ Build Vandermode vector of size @p k for element @p elt """
    result = Vector(n=k)
    result[0] = 1
    acc = 1
    for i in range(1, k):
        acc = gf_mul(acc, elt, red_pol, elt_size, elt_size, red_pol_size)
        result[i] = acc
    return result

# generator polynomial is X^8++X¨4+X¨3+X^2+1 = 0x11d (X^8 + 0d'45')
class ErasureCoder:
    def __self__(k, m):
        assert m == 2
        # number of data devices
        self.k = k
        # number of coding devices
        self.m = m
        self.CST_VECTOR_PARITY = [
            Vector([0] * k),
            build_gen_vector(k, 2)
        ]


    def encode(self, vector):
        assert vector.n == self.k
        return Vector([gf256_dotprod(vector, CST_VECTOR_PARITY[i]) for i in range(self.m)])

    def encode_alternate(self, vector):
        # build bit matrix
        raise NotImplementedError

def gf256_dotprod(v0, v1):
    assert v0.n == v1.n
    RED_POL = 0x1d
    acc = 0
    for i in range(v0.n):
        acc ^= gf_mul(v0[i], v1[i], RED_POL, 8, 8, 8)
    return acc

class Vector:
    """ Generic vector structure """
    def __init__(self, data=None, n=None):
        """ Vector init from data list or size @p n
            if only @p n is set, the vector is filled with 0s """
        if data:
            # data copy
            self.data = [v for v in data]
            self.n = len(self.data)
        elif n:
            self.data = [0] * n
            self.n = n
        else:
            raise NotImplementedError

    def __getitem__(self, index):
        """ Extract an item at index @p index """
        return self.data[index]

    def __add__(self, b):
        """ element-wise addition of 2 vectors """
        assert self.n == b.n
        return Vector([x ^ y for x, y in zip(self.data, b.data)])

    @property
    def value(self):
        """ return the vector aggregated value """
        acc = 0
        for i in range(self.n):
            acc ^= self.data[i] << i
        return acc

    def __mul__(self, v):
        """ Vector by Matrix multiplication """
        if isinstance(v, BitMatrix):
            return (BitMatrix(self.value, 1, self.n) * v).vector_row(0)
        else:
            raise NotImplementedError


TEST_NUM = 100
RED_POL = 0x1d # X^8+x^4+X^3+x^2+1
for t in range(TEST_NUM):
    a = random.randrange(2**8)
    b = random.randrange(2**8)
    meth0_result = gf_mul(a, b, RED_POL, 8, 8, 8)
    meth1_result = gf_mul_by_cst(a, b, RED_POL, 8, 8, 8)
    print(meth0_result, meth1_result)
    assert meth0_result == meth1_result
