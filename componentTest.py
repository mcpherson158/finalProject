import numpy as np

def magnitude(vector):
    x, y = vector
    magnitude = np.sqrt(x**2 + y**2)
    return magnitude

def collisionCalculaiton():
    v1 = np.array((1,0))
    v2 = np.array((0,0))
    m1 = 1
    m2 = 1
    x1 = np.array((.5,0))
    x2 = np.array((.6,0))

    massTerm1 = (2*m2)/(m1+m2)
    massTerm2 = (2*m1)/(m1+m2)
    vsub1 = v1-v2
    vsub2 = v2-v1
    xsub1 = x1-x2
    xsub2 = x2-x1
    xmag1 = magnitude(xsub1)
    xmag2 = magnitude(xsub2)
    xmagsqr1 = xmag1**2
    xmagsqr2 = xmag2**2
    dotprod1 = np.dot(vsub1,xsub1)
    dotprod2 = np.dot(vsub2,xsub2)

    prodTerm1 = dotprod1/xmagsqr1
    prodTerm2 = dotprod2/xmagsqr2


    # v1 = v1 - massTerm1*prodTerm1*xsub1
    # v2 = v2 - massTerm2*prodTerm2*xsub2

    v1 = v1 - massTerm1*prodTerm1*xsub1
    #p1.v =  v1 +v2


    #p2.v = - v2
    v2 = v2 - massTerm2*prodTerm2*xsub2

    print(v1,v2)

collisionCalculaiton()
