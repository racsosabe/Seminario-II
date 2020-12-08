from mpmath import *
import sys

def norm2(a):
    return a.real**2 + a.imag**2

def getMin(a, b):
    prod = a * b.conjugate()
    if prod.real >= norm2(a): return abs(a)
    if prod.real >= norm2(b): return abs(b)
    return abs(prod.imag / abs(b-a))

def f(x, n):
    ans = 1
    for i in range(2,n+1):
        ans += exp(-x * log(i))
    return ans

def M(s1, s2, n):
    ans = 0
    for i in range(2,n+1):
        ans += exp((-min(s1.real,s2.real)) * log(i)) * log(i) * log(i)
    return ans

def isGreater(s1, s2, localM, n):
    L = getMin(f(s1, n),f(s2, n))
    R = localM * norm2(s1 - s2) / 8
    EPS = mpf('10') ** (-20)
    return L > R + EPS

def varArg(s1, s2, n):
    EPS1 = mpf('10') ** (-25)
    if abs(s1.real - 2) <= EPS1 and abs(s2.real - 2) <= EPS1:
        return arg(f(s2, n)) - arg(f(s1, n))
    localM = M(s1,s2,n)
    t0 = mpf(0)
    t = mpf(1)
    L = s1
    d = s2 - s1
    ans = mpf(0)
    EPS2 = mpf('10') ** (-80)
    EPS3 = mpf('10') ** (-40)
    while t0 + EPS1 < 1:
        R = L + t * d
        while (not isGreater(L,R,localM,n)) and t > EPS2:
            t /= 2
            R = L + t * d
        if t <= EPS2:
            if s1.real == s2.real:
                d = EPS3 if s1.imag < s2.imag else -EPS3
                return varArg(s1 + mpc(d,0), s2 + mpc(d,0), n)
            else:
                d = EPS3 if s1.real > s2.real else -EPS3
                return varArg(s1 + mpc(0,d), s2 + mpc(0,d), n)
        ans += arg(f(R, n) / f(L, n))
        t0 += t
        L += d * t
        t = min(2 * t, 1 + EPS1 - t0)
    return ans

def getRectangle(LD, RU, n):
    LU = mpc(LD.real,RU.imag)
    RD = mpc(RU.real,LD.imag)
    return (varArg(LD,RD,n) + varArg(RD,RU,n) + varArg(RU,LU,n) + varArg(LU,LD,n)) / (2 * acos(-1))


mp.dps = 100
n = int(sys.argv[1])
values = input().split()
LD = mpc(mpf(values[0]), mpf(values[1]))
values = input().split()
RU = mpc(mpf(values[0]), mpf(values[1]))
print(getRectangle(LD, RU, n))
