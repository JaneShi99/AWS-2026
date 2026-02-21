R.<x> = PolynomialRing(ZZ)

p = 101

f = -(x^7 - x^6 + 6*x^5 - 7*x^4 + 5*x^3 + x^2 - x + 1)

C = HyperellipticCurve(f)

def count_points():
    for k in range(1, 4):
        q = p^k
        Fq = GF(q)
        Rq.<X> = PolynomialRing(Fq)
        f_q = f.change_ring(Fq)
        C = HyperellipticCurve(f_q)
        print(f"#C(F_{{{q}}}) = {C.count_points(1)[0]}")

def genus3Zeta(N1, N2, N3):
    powerPoly = N1*x + N2*x^2/2 + N3*x^3/3
    zetaFunc = 1 + (powerPoly)^1 + (1/(factorial(2)))*(powerPoly)^2 + (1/(factorial(3)))*(powerPoly)^3

    c1 = zetaFunc[1]
    c2 = zetaFunc[2]
    c3 = zetaFunc[3]
     
    Lpoly = zetaFunc*(1-x)*(1-p*x)
    print(Lpoly)

    '''
    We ended up getting: 15163*x^3 + 1167*x^2 + 51*x + 1
    '''

    return Lpoly

#ans = genus3Zeta(153, 9935, 1029891)

a1, a2, a3 = 51, 1167, 15163


PS.<T> = PowerSeriesRing(QQ, default_prec=15)

Lpoly = p^3*T^6 + p^2*a1*T^5 + p*a2*T^4 + a3*T^3 + a2*T^2 + a1*T + 1 


Zeta_exp2 = log(Lpoly/((1-T)*(1-p*T)))
print(Zeta_exp2)

'''
Answer: 153*T + 9935/2*T^2 + 343297*T^3 + 104091107/4*T^4 + ...
'''
