from galois_field import GFpn

# Generating the field GF(2^4)
# irreducible polynomial. (in this case, x^4 + x+1)
import sys

a = [1, 1, 0]
b = [0, 1, 0]
p = [
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 1, 1, 1
]

# if (len(sys.argv) > 1):
#     a = eval("[" + sys.argv[1].replace(" ", "") + "]")
# if (len(sys.argv) > 2):
#     b = eval("[" + sys.argv[2].replace(" ", "") + "]")
# if (len(sys.argv) > 3):
#     p = eval("[" + sys.argv[3].replace(" ", "") + "]")
try:
    gf = GFpn(2, p)
except Exception as e:
    print("Error:", e)
    sys.exit()

# Generating an element in GF(2^4)
aval = gf.elm(a)  # x^2+x+1
bval = gf.elm(b)  # x

# Arithmetic operations
# aval + bval  # 1x^2 + 1
# aval - bval  # 1x^2 + 1
# aval * bval  # 1x^3 + 1x^2 + 1x
# aval / bval  # 1x^3 + 1x

print("a:\t\t", aval)
print("b:\t\t", bval)
# print("b^{-1}: ", bval.inverse())
print("p:\t", gf, p)

print("\nAdd:\t\t", aval + bval)
# print("Subtract:\t", aval - bval)
print("Multiply:\t", aval * bval)

print("Power:\t ", aval ** 50)
# print("Divide:\t\t", aval / bval)
