R.<x> = PolynomialRing(ZZ)

#CM matrix method

#p = 101
#p = 184753
#p = 101
p = 6299
#f = -(x^8 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)
f = -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)
h = f^((p-1)/2)
g = (f.degree()-1)//2
Af = Matrix(QQ, g, g, lambda i, j: h[(i+1)*p-(j+1)])
print(Af.apply_map(lambda x: Integer(x) % p))

for r in range(1, 4):
    print("#C(F_(p^r)) mod p")
    print((1-(Af^r).trace()) % p)

'''
print(153%101, 9935%101, 1029891%101)

#Very slow method
print("\n#Very slow method")
for r in [1,2]: #ignore r=3
    hr = f^((p^r-1)/2)
    solnsModP = 1-sum([hr[i*(p^r-1)] for i in range(1,g+1)])
    print(r)
    print(solnsModP % p)
'''