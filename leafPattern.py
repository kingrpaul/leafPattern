"""
author: king.r.paul@gmail.com
"""
def abutted(a, b):
    TOLERANCE = 0.1
    if abs(a + b) < TOLERANCE: 
        return True
    else: 
        return False

def effSq(seg):
    """Returns a weighted average effective field size from the leaf pattern, 
    'seg', where 'seg' is a list of tuples (A-leaf, B-leaf) whose elements 
    give the distance of the leaf from from the center-line. """ 

    MILLENIUM, BRAINLAB = 120, 52 
    numLeaves = 2*len(seg)

    # T: leaf thicknesses by leaf pair for MLC models
    if   numLeaves == MILLENIUM:
        T = [1.0 for i in range(10)] + \
            [0.5 for i in range(40)] + \
            [1.0 for i in range(10)]
    elif numLeaves == BRAINLAB:
        T = [0.55 for i in range(3 )] + \
            [0.45 for i in range(3 )] + \
            [0.3  for i in range(14)] + \
            [0.45 for i in range(3 )] + \
            [0.55 for i in range(3 )]

    # Y: y component of distance from (0,0) to leaf center by leaf pair
    Y = [(-0.5*sum(T) + sum(T[:i]) + 0.5*T[i])  for i in range(len(T))]

    # effective field size from effective width and legnth
    area, effX, effY, numer, denom = 0.0, 0.0, 0.0, 0.0, 0.0
    for i in range(len(seg)):
        A, B = seg[i][0], seg[i][1]
        # zero for closed leaf-pairs
        if not abutted(A, B):
            area += T[i]*(A+B)
            # zero for leaf past mid-line
            A, B = max(0.0, A), max(0.0, B)
            distSqrA = Y[i]**2 + A**2
            distSqrB = Y[i]**2 + B**2
            numer +=  T[i] * A/distSqrA +  T[i] * B/distSqrB
            denom += (T[i]    /distSqrA) +(T[i]    /distSqrB)
    try:
        effX = 2.0 * numer / denom  
        effY = area / effX
        return 2.0*effX*effY/(effX+effY)
    except ZeroDivisionError: return 0.0 
    
def test():
    shape = [
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (3.03, 2.47),
        (2.88, 2.46),
        (3.08, 2.51),
        (2.86, 2.46),
        (2.88, 2.46),
        (2.91, 5.04),
        (2.5, 5.04),
        (2.55, 4.87),
        (2.38, 4.61),
        (2.38, 7.04),
        (2.61, 7.46),
        (2.48, 6.55),
        (3.02, 6.52),
        (3.9, 7.2),
        (4.5, 7.5),
        (4.5, 7.5),
        (4.5, 7.5),
        (4.5, 7.5),
        (4.45, 7.5),
        (4.0, 7.5),
        (3.5, 7.5),
        (3.49, 7.5),
        (3.0, 7.5),
        (3.0, 7.5),
        (3.0, 7.5),
        (2.5, 7.5),
        (2.5, 7.5),
        (2.49, 6.52),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, 0.0)]

    print 'effSq Test Passed: ', abs(effSq(shape) - 10.725) < 0.05
    
if __name__ == '__main__': 
    test()
    
