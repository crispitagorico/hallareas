import numpy as np
import itertools, functools, sys
from free_lie_algebra import *

def labelTuple2String(lab, symbol='area'):
    assert type(lab) is tuple
    if len(lab)==1:
        return str(lab[0])
    return symbol+"("+labelTuple2String(lab[0],symbol)+","+labelTuple2String(lab[1],symbol)+")"

def label2String(lab,symbol='area'):
    if len(lab)==1:
        return labelTuple2String(lab[0],symbol)
    return " "+(" * ".join(labelTuple2String(i,symbol) for i  in lab))

def PolynomialRegression(element, monomials, symbol, width, T):
    
    label = element.pretty()
    letters_ = tuple(range(1, width+1))
    letters = tuple(letter2Elt(i) for i in letters_)

    # Polynomial regression
    coeffs = T.inTermsOf(element, [j for i,j in monomials])
    labels = [k[0] for k in monomials]
    
    basis_elements = []
    foliages = {}
    for lab, c in zip(labels, coeffs):
        if 1e-10<math.fabs(c):
            basis_elements.append('[{:.3f}] {}'.format(c, label2String(lab,symbol)))
            for l in lab:
                foliages[labelTuple2String(l,symbol)] = c
    summation = ''
    for k in basis_elements:
        summation += '  ' + k + '  +'
    
    return summation[:-1]

# given a tuple of pairs, return list of the firsts and shuffle of seconds
def shufflePairs(a):
    return list(i[0] for i in a),functools.reduce(shuffleProduct,(i[1] for i in a))

def lieFromTree(tup):
    assert type(tup) is tuple
    if len(tup)==1:
        return(letter2Elt(tup[0]))
    assert len(tup)==2
    return lieProduct(lieFromTree(tup[0]),lieFromTree(tup[1]))

def hsFromTree(tup):
    assert type(tup) is tuple
    if len(tup)==1:
        return(letter2Elt(tup[0]))
    assert len(tup)==2
#     return rightHalfShuffleProduct(hsFromTree(tup[0]),hsFromTree(tup[1]))
    return leftHalfShuffleProduct(hsFromTree(tup[0]),hsFromTree(tup[1]))


def area(a,b):
#     return rightHalfShuffleProduct(a,b)-rightHalfShuffleProduct(b,a)
    return leftHalfShuffleProduct(a,b)-leftHalfShuffleProduct(b,a)

def area_(a):
    return functools.reduce(area,a)

def areaFromTree(tup):
    assert type(tup) is tuple
    if len(tup)==1:
        return(letter2Elt(tup[0]))
    assert len(tup)==2
    return area(areaFromTree(tup[0]),areaFromTree(tup[1]))


#area_ of all the combinations of depth letters 
def leftBracketedAreasFromLetters(letters,depth):
    assert depth>1
    return [area_(letter2Elt(j) for j in i) for i in itertools.product(letters,repeat=depth) if i[0]<i[1]]

def leftBracketedAreasFromDistinctLetters(letters):
    return [area_(letter2Elt(j) for j in i) for i in itertools.permutations(letters) if i[0]<i[1]]

def allAreasFromLetters(letters,depth):
    gradedlist = [letters]
    for newDepth in range(2,depth+1):
        gradedlist = gradedlist + [[]]
        for leftDepth in range(1,newDepth):
            rightDepth=newDepth-leftDepth
            for left in gradedlist[leftDepth-1]:
                for right in gradedlist[rightDepth-1]:
                    gradedlist[-1].append(area(left,right))
    return gradedlist[-1]

def conjecturedBasis(letters,depth):
    assert depth>1
    out=[]
    lets=[letter2Elt(i) for i in letters]
    for final2 in itertools.combinations(lets,2):
        f=lieProduct(*final2)
        for rest in itertools.product(lets,repeat=depth-2):
            out.append(concatenationProductMany(rest+(f,)))
    return out