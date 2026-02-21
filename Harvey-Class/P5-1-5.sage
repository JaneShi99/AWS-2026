def random_poly(g, p):
    R.<x> = PolynomialRing(GF(p))
    while True:
        coeffs = [GF(p).random_element() for _ in range(2*g + 2)]
        f = R(coeffs + [1])
        if f.discriminant() != 0:
            return f
    
    
def get_h_k(f, k, p, h):
    d = f.degree()
    return (1/(k*f[0])) * sum(((GF(p)(j)/2) - k) * f[j] * h[k-j] for j in range(1, d))

def compute_h_recursively(f, k, p):
    R.<x> = PolynomialRing(GF(p))
    if k == 0:
        return R(f[0]^((p-1)/2))
    else:
        h_prev = compute_h_recursively(f, k-1, p)
        return h_prev + get_h_k(f, k, p, (h_prev))*x^k


def experiment():
    p = 11
    g = 5 
    f = random_poly(g, p)

    R.<x> = PolynomialRing(GF(p))

    h_computed = compute_h_recursively(f, p-1, p)
    print(h_computed)
    print("\n")

    h_poly = f^((p-1)/2)
    print(h_poly)

experiment()
