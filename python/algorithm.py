import numpy as np
import multiprocessing
from multiprocessing import Pool
import itertools
from os import listdir
from os.path import isfile, join
import sympy
import solmodn

#fixed points are represented by [(a_1, ..., a_dim), d], all integers with d non-zero, which corresponds to (exp(i*a_1/d), ..., exp(i*a_dim/d))

#given a matrix g and a point pd checks if (g + I) fixes pd
def fix(g, pd):
    p = pd[0]
    d = pd[1]    
    for n in g.dot(np.array(p)):
        if n % d != 0:
            return 0
    return 1
    
'''
Given an invertible matrix A, returns fixed points of A + I
By considering action of g^{-1} on the points of Z^n we can see that each fixed point has to have the form [(a_1, ..., a_dim), det(matrix)].
'''
def fixedpoints(A):    
    d = int(round(np.abs(np.linalg.det(A))))    
    return [[x,d] for x in solmodn.solutions(A, d)]

'''
Checks if the group is excluded based on Lemma 4.3 with k=1 and k=2.
To check if a point is fixed by all matrices in the group, it is enough to check if all generators of the group fix it.
If a group is excluded based on these conditions, it creates a .sol file with the value of k, a matrix or a pair of matrices(depending on the value of k) which have a trivial intersection of kernels,
a point which is fixed by them, and one of the generators of the group that doesn't fix it.
'''
def check_group(task):
    
    ff = task[0]    
    group = task[1]
    gens = task[2]
    dim = task[3]
    
    #k=1
    det1 = [g for g in group if int(round(np.linalg.det(g - np.identity(dim, dtype = np.intc)))) != 0]    
    for g in det1:
        for point in fixedpoints(g - np.identity(dim, dtype = np.intc)):
            for gen in gens:
                if fix(gen - np.identity(dim, dtype = np.intc), point) == 0:
                    outfile = open("./dim{}/".format(dim)+ff[:-4]+".sol", "w")
                    outfile.write("k={}\ng1={}\ngen={}\npoint={}".format(1, g, gen, point))
                    outfile.close()  
                    return [ff, 1]
                
    #k=2 
    i=0
    pairs = itertools.combinations(group[2:], 2)    
    for (g1,g2) in pairs:        
            '''We check if ker(g-I) \cap ker(g2 - I) is empty by building a matrix consisting of rows of g and g2, and computing its rank. If it has full rank, we extract a set of linearly independent rows into the matrix comb. Then common fixed points of g, g2 are certainly also fixed points of comb + I'''
            comb = np.concatenate((g1 - np.identity(dim, dtype = np.intc),g2 - np.identity(dim, dtype = np.intc)), dtype = np.intc)
            _, inds = sympy.Matrix(comb).T.rref()
            
            if len(inds) == dim:                
                sliced = comb[list(inds)]
                fixed = fixedpoints(sliced)
                fixed = [point for  point in fixed if fix(g1 - np.identity(dim, dtype = np.intc), point) == 1 and fix(g2 - np.identity(dim, dtype = np.intc), point) == 1]                
                for point in fixed:
                    for gen in gens:
                        if fix(gen - np.identity(dim, dtype = np.intc), point) == 0:
                            outfile = open("./dim{}/".format(dim)+ff[:-4]+".sol", "w")                           
                            outfile.write("k={}\ng1={}\ng2={}\ngen={}\npoint={}".format(2, g1, g2, gen, point))
                            outfile.close()        
                            return [ff, 1, point]           

    return [ff,0]

#converts a string representation of group elements to a list of matrices
def stringtogroup(s, dim):    
    s = s.translate({ord('['): None, ord(']'): None, ord('{'): None, ord('}'): None, ord(','): " ", ord('\\'): None, ord('\n'): None})
    arr = [int(x) for x in s.split()]    
    m = [arr[i:i + dim * dim] for i in range(0, len(arr), dim * dim)]
    return([np.array([x[i:i + dim] for i in range(0, len(x), dim)], dtype = np.intc) for x in m])

#given a filename and dimension, returns a list of matrices
def filetogroup(file, dim):
    with open (file, "r") as myfile:
        data = myfile.read()
        return stringtogroup(data, dim)


if __name__ == "__main__":
    for dim in [4, 5, 6]:
        print("Dimension: {}".format(dim))
        
        path = "./dim{}/".format(dim)
        files = listdir(path)        
        groups = [f for f in files if f.endswith(".grp")]
               
        tasks = []                
        for f in groups:        
            if (f[:-4] + ".sol") not in files:
                genspath = path + f[:-4] + ".gens"
                gens = filetogroup(genspath, dim)
                grp = filetogroup(path + f, dim)
                tasks.append([f, grp, gens, dim])                
                
        print("{} groups to test".format(len(tasks)))       
        
        with multiprocessing.Pool() as pool:
            done = []
            for res in pool.imap_unordered(check_group, tasks):                
                if res[1] != 0:                    
                    done.append(res[0])                
        print("{} groups excluded".format(len(done)))
