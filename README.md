# gftk
Galois Field ToolKit

set of python utility to implement computation on Galois Field of characteristic 2.
This is intended to compare BitMatrix based methods with standard methods.

# HowTo

Tested with python3 (3.5.3)

```
import utility
```


## Carry Less Multiplication

To compute a multiplication between two polynomials of GF(2){X]
simply use the `clm` function with the polynomials as parameter (encoded in binary representation).
No extra parameters is required.

```python
from utility import clm

# random left hand side 8-bit input
a = 0xfe

# random right hand side 8-bit input
b = 0x17

c = clm(a, b)

print("{:x}".format(c))
```
will produce the following output:
```
0xd1a
```

## Galois Field Multiplication

The funcion `gf_mul` actually implements a polynomial multiplication
followed by a modulo reduction.
Multiply operands and modulo polynomial must be given as parameters along their respective sizes.
The modulo polynomial is expressed without the most significant monomyals.
For example `0x1d` with associated size 8 encodes X^8+X^4+X^3+X^2+1.


```python
from utility import gf_mul

# random left hand side 8-bit input
a = 0xfe
a_size = 8 # a's size in bits

# random right hand side 8-bit input
b = 0x17
b_size = 8 # b's size in bits

# polynomial use for reduction (without most significant monomials)
RED_POL = 0x1d # (X^8)+X^4+X^3+X^2+1
mod_size = 8 # reduction polynomial is 8-bit wide (omitting MSB)

c = gf_mul(a, b, RED_POL, a_size, b_size, mod_size)

print("{:x}".format(c))
```
will produce the following output:
```
0x9b
```


## Bit Matrices

```python
from utility import BitMatrix

# build a bit-matrix 8x8 full of zeros
zeros = BitMatrix(0, row=8, col=8)

# build a 8x8 identity matrix
identity = BitMatrix.identity(8)

# transpose the identity matrix (giving identity ...)
id_transposed = identity.transpose()

# A numerical value can be used to initalize the bit matrix elements
# the value is assumed to encode in a row major and little endian layout
# the matrix

# for example for a 8x8 matrix will the first row full of '1':
first_full_row = BitMatrix(0xff, row=8, col=8)
print("first_full_row:\n{}".format(str(first_full_row)))

# for example for a 8x8 matrix will the first column full of '1':
# each 0x01 encodes a '1' in the first column of a different rows
# the leftmost 0x01 is the 7-th row (because it is the 7-th byte of
# the numerical constant)
first_full_col = BitMatrix(0x0101010101010101, row=8, col=8)
print("first_full_col:\n{}".format(str(first_full_col)))
```

```
first_full_row:
1 1 1 1 1 1 1 1
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0

first_full_col:
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0
```
