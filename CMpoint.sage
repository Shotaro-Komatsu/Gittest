def period_matrix(phis, ideal):
    basis = ideal.basis()
    riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
    symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))

    big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
    big_period_matrix.subdivide(2,2)

    Z1 = big_period_matrix.subdivision(0,0)
    Z2 = big_period_matrix.subdivision(0,1)
    Z  = Z2.inverse()*Z1
    return Z
