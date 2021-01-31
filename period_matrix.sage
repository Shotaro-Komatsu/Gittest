def period_matrix(bas, embs):
    M = Matrix([[phi(b) for b in bas] for phi in embs])
    Omega1 = M[:,:2]
    Omega2 = M[:,2:]
    return Omega2.inverse()*Omega1

def split(bas):
    g = ZZ(bas.degree()/2)
    return vector(bas[:g]), vector(bas[g:])

def join(bas1,bas2):
    return vector(list(bas1)+list(bas2))

def B_to_M(B):
    g = B.nrows()
    assert B.ncols() == g
    return ABCD_to_M(identity_matrix(g), B, zero_matrix(g, g), identity_matrix(g))

def ABCD_to_M(A,B,C,D):
    g = A.ncols()
    for E in [A,B,C,D]:
        assert E.nrows() == E.ncols() == g
    return Matrix([join(A[i],B[i]) for i in range(A.nrows())] +
                  [join(C[i],D[i]) for i in range(A.nrows())])
