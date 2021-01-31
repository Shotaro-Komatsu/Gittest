load('/home/sage/Documents/ドライブ/Sage/code/CMfield.sage')
load('/home/sage/Documents/ドライブ/Sage/code/theta_function.sage')
load('/home/sage/Documents/ドライブ/Sage/code/CMpoint.sage')

K=CMFieldfromPoly(2,-5,'a',prec=200)

O=K.suitable_ideals()[1]
conj=K.complex_conjugation()
principal=K.different()*O*(conj(O))


pair=K.xi_and_CMtype()

delta = pair[0][0]
phis = [pair[0][1][1], pair[0][1][0]]
xi=delta^(-1)

Omega = period_matrix(phis,O)

jlist = Igusa_invariants(Omega,200)
