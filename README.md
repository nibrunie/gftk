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
