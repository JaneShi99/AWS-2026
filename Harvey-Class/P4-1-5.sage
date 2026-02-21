R.<x> = PolynomialRing(ZZ)
p = 101
f = -(x^7 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)

h = f^((p-1)/2)
g = 3
solnsModP = 1 - sum([h[i*(p-1)] for i in range(1,g+1)])
print(solnsModP % 101)
 
print(153 % 101)