# -*- coding: utf-8 -*-


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

def bit_matrix_multiply(a, b):
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
        return self.matrix[row_id]

    def __mul__(self, b):
        return bit_matrix_multiply(self, b)

    @staticmethod
    def identity(row_col):
        result = BitMatrix(0, row_col, row_col)
        for i in range(row_col):
            result[i][i] = 1
        return result

    @property
    def transpose(self):
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

    def __str__(self):
        return "\n".join(" ".join("{}".format(self.matrix[i][j]) for j in range(self.col)) for i in range(self.row))



