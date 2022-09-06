import numpy as np
import itertools
from functools import reduce

'''
For a d x d matrix A, a prime p and a natural number k, finds all solutions x \in (Z/p^kZ)^d to equation Ax=0
'''
def solutionsprimepower(A, p, k):
    d = np.shape(A)[0]    
    a = [range(0, p)] * d
    t = [np.array(x) for x in list(itertools.product(*a))]
    if k == 1:
        sols=[]
        for y in t:            
            r = A@y
            sol = True
            for n in r:
                if n%p != 0:
                    sol = False
            if sol is True:
                sols.append(y)
        return sols
    else:
        prevsols = solutionsprimepower(A, p, k-1)        
        sols=[]
        for prevsol in prevsols:            
            tt = [prevsol + p**(k-1) * z for z in t]            
            for y in tt:
                r = A@y
                sol = True
                for n in r:
                    if n%(p**k) != 0:
                        sol = False
                if sol is True:
                    sols.append(y)
        return sols

#The next two functions were taken from https://rosettacode.org/wiki/Chinese_remainder_theorem#Python
#Computes the modular inverse of a mod b
def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while a > 1:
        q = a // b
        a, b = b, a%b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1

#For n_1, ..., n_k pairwise coprime, returns the unique solution in [0, n_1* ...* n_k] of a system of equations x_i = a_i mod n_i, i = 1, ... ,k. 
def chinese_remainder(n, a):
    sum = 0
    prod = reduce(lambda a, b: a*b, n)
    for n_i, a_i in zip(n, a):
        p = prod // n_i
        sum += a_i * mul_inv(p, n_i) * p
    return sum % prod 


#returns factorization of a natural number n in the form [[p1, a1], [p2, a2], ...]
def factorize(n):
    k = 2
    factors = []
    while n != 1:
        count = 0
        while n%k == 0:
            n/=k
            count+=1            
        if count != 0:
            factors.append([k, count])
        k+=1
    return factors

#For a d x d matrix A, and a natural number n, finds all solutions x \in (Z/nZ)^d to equation Ax=0
def solutions(A, n):
    if n==1:
        return np.zeros(np.shape(A))
    sols = []
    d = np.shape(A)[0]
    factors = factorize(n)
    factorsols = [solutionsprimepower(A, f[0], f[1]) for f in factors]
    combs = list(itertools.product(*factorsols))
    for comb in combs:
        n = [a[0] ** a[1] for a in factors]
        sol = np.zeros([d])
        for i in range(0, d):
            a = [x[i] for x in comb]
            sol[i] = chinese_remainder(n, a)
        sols.append(sol)
    return sols

    
