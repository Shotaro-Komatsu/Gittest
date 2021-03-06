load('CMfield.sage')
load('theta_function.sage')
load('CMpoint.sage')

K=CMFieldfromPoly(2,-5,'a',prec=200)

a=K.suitable_ideals()[1]
conj=K.complex_conjugation()
principal=K.different()*a*(conj(a))


pair=K.xi_and_CMtype()

delta = pair[0][0]
phis = [pair[0][1][1], pair[0][1][0]]
xi=delta^(-1)

Omega = period_matrix(phis,O)

jlist = Igusa_invariants(Omega,200)
