import free_lie_algebra as fla

def larea(a,b): #from Joscha's conjectures.conjecture_1()
    return fla.leftHalfShuffleProduct(a,b) - fla.leftHalfShuffleProduct(b,a)
    
def lvol(a,b,c):
    return larea(larea(a,b),c)+larea(larea(b,c),a)+larea(larea(c,a),b)
    
shuffleProduct=fla.shuffleProduct
lhs=fla.leftHalfShuffleProduct

x=fla.parse('1')
y=fla.parse('2')
z=fla.parse('3')
w=fla.parse('4')

print('Shuffle pullout')
print((3*larea(z,shuffleProduct(x,y))).pretty())
print((shuffleProduct(x,larea(z,y))+shuffleProduct(y,larea(z,x))-shuffleProduct(x,shuffleProduct(y,z))+larea(larea(z,y),x)+larea(larea(z,x),y)).pretty())

print('Area-Jacobi')
print((larea(larea(x,y),z)+larea(larea(y,z),x)+larea(larea(z,x),y)).pretty())
print((shuffleProduct(x,larea(y,z))+shuffleProduct(y,larea(z,x))+shuffleProduct(z,larea(x,y))).pretty())

print('Tortkara 1')
print(larea(larea(x,y),larea(x,z)).pretty())
print(larea(x,lvol(x,y,z)).pretty())

print('Tortkara 2')
print((larea(larea(x,y),larea(w,z))+larea(larea(z,y),larea(w,x))).pretty())
print((larea(x,lvol(y,z,w))+larea(z,lvol(y,x,w))).pretty())

print((2*lhs(lhs(z,lhs(z,y)),lhs(x,lhs(z,lhs(z,x))))+lhs(lhs(z,x),lhs(y,lhs(z,lhs(z,lhs(z,x))))) +lhs(lhs(z,x),lhs(lhs(z,x),lhs(z,lhs(z,y)))) -lhs(lhs(z,y),lhs(x,lhs(z,lhs(z,lhs(z,x))))) -lhs(lhs(z,y),lhs(lhs(z,x),lhs(z,lhs(z,x)))) -lhs(lhs(z,lhs(z,x)),lhs(x,lhs(z,lhs(z,y)))) -lhs(lhs(z,lhs(z,x)),lhs(y,lhs(z,lhs(z,x)))) -lhs(lhs(z,lhs(z,lhs(z,x))),lhs(x,lhs(z,y))) +lhs(lhs(z,lhs(z,lhs(z,x))),lhs(y,lhs(z,x)))))

