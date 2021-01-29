load('/home/sage/Documents/ドライブ/Sage/code/CMFIELD.sage')
load('/home/sage/Documents/ドライブ/Sage/code/theta_funct.sage')
D= 6
P = Primes()
Num=0
Ilists=[]

E = -3
K=CMFieldfromPoly(D,E,'a',prec=200)
K0=K._K0
h0=K0.class_number()
h=K.class_number()
rho,tau=K0.embeddings(K)
L=K.relativize(rho(K0.gen()),'a')
L_into_K, K_into_L = L.structure()
L1=L.base_field()
OL=L.maximal_order()
a0=L.gen()
a1=L1.gen()

O=K.suitable_ideals()[Num]
conj=K.complex_conjugation()
principal=K.different()*O*(conj(O))


allphis=K.all_CMtypes()
pair=[]
for Phis in allphis:
    delta=K.good_generator(Phis,principal)
    if not (delta == None):
        pair.append([delta, Phis])

delta = pair[0][0]
phis = pair[0][1]

xi=delta^(-1)
basis = O.basis()

riemann_form = Matrix(ZZ,[[(conjugate(x)*xi*y).trace() for y in basis] for x in basis])
symplectic_basis = Sequence(riemann_form.symplectic_form()[1]*vector(basis))
big_period_matrix = Matrix([[phi(b) for b in symplectic_basis] for phi in phis])
big_period_matrix.subdivide(2,2)
Omega1 = big_period_matrix.subdivision(0,0)
Omega2 = big_period_matrix.subdivision(0,1)

Omega=Omega2.inverse()*Omega1

Igusalist=Igusa_invariants(Omega,50)


'''
intIlist=[round(Igusalist[0].real()), round(Igusalist[1].real()), round(Igusalist[2].real())]
intfIlist=[factor(intIlist[0]), factor(intIlist[1]), factor(intIlist[2])]
'''


'''
abas=F.gen(0)
a0bas=-1/2*K0.gen(0)
basphi=K.relative_basis()
Phi=basphi[0][0]
tau=basphi[1][0]
xi=K.xi_by_basis()


Omega1=Matrix([[Phi[0](tau*omega),Phi[0](tau)],[Phi[1](tau*omega),Phi[1](tau)]])
Omega2=Matrix([[Phi[0](1),Phi[0](omega)],[Phi[1](1),Phi[1](omega)]])
Omega=Omega2.inverse()*Omega1
'''
