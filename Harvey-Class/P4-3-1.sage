R.<x> = PolynomialRing(ZZ)

def countModp(f, p, r):
    h = f^((p-1)/2)
    g = (f.degree()+1)//2
    Af = Matrix(QQ, g, g, lambda i, j: h[(i+1)*p-(j+1)])
    return (1-(Af^r).trace()) % p

def countPoints(f, p, r):
    Fq = GF(p^r)
    fFq = f.change_ring(Fq)
    # Use quadratic character: 1 if square, -1 if non-square, 0 if zero
    def quad_char(val):
        if val == 0:
            return 0
        return 1 if val.is_square() else -1
    return sum(1 + quad_char(fFq(a)) for a in Fq) + 1

f1 = x^7 + 1
f2 = x^7 - x^5 + 1
p = 3

for r in range(1, 4):
    print("f1")
    print(countModp(f1, p, r))
    print("f2")
    print(countModp(f2, p, r))
    
def zetaFunc(f, p, prec):
    PS.<T> = PowerSeriesRing(QQ, default_prec=prec)
    Zeta = exp(sum(countPoints(f, p, r)*T^r/r for r in range(1, prec+1)))
    Lpoly = Zeta*(1-T)*(1-p*T)
    return Lpoly

print(zetaFunc(f1, 3, 4))
print(zetaFunc(f2, 3, 4))


def getAf(f, p):
    h = f^((p-1)/2)
    g = (f.degree()+1)//2
    Af = Matrix(QQ, g, g, lambda i, j: h[(i+1)*p-(j+1)])
    return Af

PS.<T> = PowerSeriesRing(QQ)

Af1 = getAf(f1, 3)
Af2 = getAf(f2, 3)

Af1T = Af1.change_ring(QQ[T])
Af2T = Af2.change_ring(QQ[T])

print((1-T*Af1T).det())
print((1-T*Af2T).det())
