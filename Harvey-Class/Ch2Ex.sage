from sage.all import getrandbits
from timeit import timeit

def ex211():
    def helper(nbits):
        a, b = getrandbits(10**nbits), getrandbits(10**nbits)
        t = timeit(lambda: a*b, number=1)
        return t
    
    for i in range(3, 8):
        print(str(i)+"-bit multiplication takes "+str(helper(i)) + " seconds")
    
    return 0

def ex213a():
    def helper(nbits):
        a, b = getrandbits(2*10**nbits), getrandbits(10**nbits)
        t = timeit(lambda: a // b, number=1)
        return t

    for i in range(3, 8):
        print(str(i)+"-bit division takes "+str(helper(i)) + " seconds")
    
def ex213b():
    def helper(nbits):
        a, b = getrandbits(10**nbits), getrandbits(10**nbits)
        t = timeit(lambda: gcd(a,b), number=1)
        return t
    
    for i in range(3, 8):
        print(str(i)+"-bit gcd takes "+str(helper(i)) + " seconds")
    
if __name__ == "__main__":
    ex211()
    print("\n")
    ex213a()
    print("\n")
    ex213b()
    