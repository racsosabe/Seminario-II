from mpmath import *
import sys
import copy

def add(a, b):
    c = []
    for i in range(len(a)):
        c.append(a[i] + b[i])
    return c

def mul(c, a):
    res = []
    for i in range(len(a)):
        res.append(c * a[i])
    return res

def dot(a, b):
    n = len(a)
    m = len(b)
    assert n == m
    ans = mpf('0')
    for i in range(n):
        ans += a[i] * b[i]
    return ans


def getGramSchmidt(v):
    n = len(v)
    bp = [mpf('0') for i in range(n)]
    for i in range(n):
        bp[i] = v[i]
        j = i - 1
        while j >= 0:
            mu = dot(v[i], bp[j]) / dot(bp[j], bp[j])
            for pos in range(len(bp[i])):
                bp[i][pos] -= mu * bp[j][pos]
            j -= 1
    return bp

def LLL(base, delta):
    b = copy.deepcopy(base)
    n = len(b)
    bp = copy.deepcopy(b)
    bp = getGramSchmidt(bp)
    mu = []
    B = []
    for i in range(n):
        B.append(dot(bp[i], bp[i]))
        j = 0
        cur = []
        while j <= i:
            cur.append(dot(b[i], bp[j]) / dot(bp[j], bp[j]))
            j += 1
        mu.append(cur)
    # Go with LLL
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q == 0: continue
            for pos in range(len(b[k])):
                b[k][pos] -= q * b[j][pos]
        
        for pos in range(k + 1):
            if pos <= k:
                mu[k][pos] = dot(b[k], bp[pos]) / dot(bp[pos], bp[pos])

        if B[k] >= (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            bp = copy.deepcopy(b)
            bp = getGramSchmidt(bp)
            for i in range(n):
                B[i] = dot(bp[i], bp[i])
                j = 0
                while j <= i:
                    mu[i][j] = dot(b[i], bp[j]) / dot(bp[j], bp[j])
                    j += 1
            if k >= 2: k -= 1
    for i in range(1, n): # Final check :v 
        assert dot(bp[i], bp[i]) >= (dot(b[i], bp[i - 1]) / dot(bp[i - 1], bp[i - 1])) ** 2 * dot(bp[i - 1], bp[i - 1])
    return b

def build(alpha, eps):
    n = len(alpha)
    b = [[mpf('0') for j in range(n + 1)] for i in range(n + 1)]
    for i in range(n):
        b[i][i] = mpf('1')
    for i in range(n):
        b[n][i] = -alpha[i]
    b[n][n] = mpf('0.5') ** (n * (n + 1) / mpf('4.0')) * eps ** (n + 1)
    return b

dim = int(sys.argv[1])
mp.dps = 100
vuelta = mpf('2') * acos(mpf('-1'))
eps1 = mpf('10') ** (-2)
eps = eps1 / vuelta
print(eps)
delta = mpf('3/4')

def isPrime(x):
    if x == 1: return False
    i = 2
    while i * i <= x:
        if x % i == 0:
             return False
        i += 1
    return True

alpha = []

for i in range(2, dim + 1):
    if isPrime(i):
        alpha.append(log(i) / vuelta)

n = len(alpha)
b = build(alpha, eps)
basis = LLL(b, delta)
p = [mpf('0') for i in range(n)]
q = mpf('2') ** (n * (n + 1) / mpf('4')) * eps ** (-n - 1) * basis[0][n]
for i in range(n):
    p[i] = basis[0][i] + q * alpha[i]

assert sqrt(dot(basis[0], basis[0])) <= eps

if q < 0:
    q = -q
    for i in range(n):
        p[i] = -p[i]
    for i in range(n + 1):
        basis[0][i] = -basis[0][i]

print(q)
