from sage.rings.number_field.number_field import NumberField_absolute
from sage.rings.number_field.number_field_rel import NumberField_relative
from sage.numerical.mip import MIPSolverException


class CMFieldfromPoly(NumberField_absolute):

    def __init__(self, D,E, names, prec = 664):
        self.polynomial = (sqrt(D)+sqrt(E)).minpoly()
        NumberField_absolute.__init__(self, self.polynomial, names)
        self._prec = prec
        for i in range(3):
            if self.subfields(2)[i][0].is_CM()==False:
                self._K0 = self.subfields(2)[i][0]
                self._embedK0toK = self.subfields(2)[i][1]
        self._CMPoints = [] #here we can cache the CM points of K

    def K0(self):
        return self._K0

    def embedK0toK(self):
        return self._embedK0toK


    def one_CMtype(self, prec=None):
        """
        KのCMタイプの１つ
        """
        if prec == None:
            prec = self._prec
        embeddings = []
        for embed in self.complex_embeddings(prec):
            embeddings.append(embed)
        phi1= embeddings.pop(0) #popは選びながら除外
        phi1bar = find_conj(embeddings,phi1)
        embeddings.pop(embeddings.index(phi1bar))
        phi2 = embeddings.pop(0)
        candidate = [phi1,phi2]
        return candidate


    def all_CMtypes(self, prec = None):
        """
        Kの全てのCMタイプのリスト
        """
        if prec == None:
            prec = self._prec
        primitives = []
        all_triples = Combinations(self.complex_embeddings(prec),2).list()
        for triple in all_triples:
            if is_CMtype(triple):
                primitives.append(triple)
        self._allCMtypes = primitives
        return self._allCMtypes


    def suitable_ideals(self):
        """
        a*bar(a)が単項になるようなaを選ぶ
        """
        ideals = [J.ideal() for J in list(self.class_group())]
        conj = self.complex_conjugation()
        suitables = []
        for ideal in ideals:
            princ_ideal = ideal * conj(ideal)
            if princ_ideal.is_principal():
                suitables.append(ideal)
        self._suitable_ideals = suitables
        return self._suitable_ideals

    def relative_field(self):
        K0=self._K0
        E =self._E
        p=sqrt(E).minpoly()
        F.<a>=K0.extension(p)
        self._relative_field = F
        return self._relative_field

    def Fone_CMtype(self, prec=None):
        """
        FのCMタイプの１つ
        """
        try:
            F = self._relative_field
        except:
            F = self.relative_field()
        if prec == None:
            prec = self._prec
        embeddings = []
        for embed in F.complex_embeddings(prec):
            embeddings.append(embed)
        phi1= embeddings.pop(0) #popは選びながら除外
        phi1bar = find_conj(embeddings,phi1)
        embeddings.pop(embeddings.index(phi1bar))
        phi2 = embeddings.pop(0)
        candidate = [phi1,phi2]
        return candidate


    def generators_U_plus(self,case = False):
        """
        U+の生成元
        """
        K0 = self._K0
        embedK0toK = self._embedK0toK
        u = UnitGroup(K0).fundamental_units()
        if K0.complex_embeddings()[0](u[0]) < 0:
            v = -u[0]
        else:
            v = u[0]
        generators = []
        case_number = ''
        if K0.complex_embeddings()[1](v) > 0:
            generators = [v]
            case_number = '1'
        else:
            generators = [v^2]
            case_number = '2'
        self._generators_U_plus = generators
        if not case:
            return self._generators_U_plus
        elif case:
            return [self._generators_U_plus,case_number]

    def in_U_1(self,u):
        """
        uがU1に含まれるかどうか判定
        """
        UK = self.unit_group()
        zeta_ord = UK.zeta_order()
        v0, v1 = UK.gens()
        conj = self.complex_conjugation()
        e10,e11 = UK.log(conj(v1))
        U_1  = span([[e10,e11+1],[zeta_ord,0]],ZZ)
        return (vector(ZZ,UK.log(u)) in U_1)

    def units_epsilon(self):
        """
        U+/U1を決定
        """
        generators, case = self.generators_U_plus(case = True)
        embedK0toK = self._embedK0toK
        if case == '2':
            return [self.one()]
        elif case == '1':
            if self.in_U_1(embedK0toK(generators[0])):
                return [self.one()]
            else:
                return [self.one(), embedK0toK(generators[0])]

    def unit_maker(self, tietze_list):
        """
        In units_for_xi we use a free group to compute a quotient of unit groups. This converts an element of the free group in Tietze form to an element of the unit group.
        """
        unit = self.one()
        for el in tietze_list:
            if el == 1:
                unit *= u0
            elif el == -1:
                unit *= u0**(-1)
            elif el == 2:
                unit *= u1
            elif el == -2:
                unit *= u1**(-1)
        return unit

    def units_for_xi(self):
        """
        Returns a list of representatives of the group U_K/U^+, where U_K are all units in K, and U^+ are all totally positive units in K0.
        """
        K0 = self._K0
        embedK0toK = self._embedK0toK
        UK = self.unit_group()
        try:
            generators_U_plus = self._generators_U_plus
        except:
            generators_U_plus = self.generators_U_plus()
        logs = [UK.log(embedK0toK(u)) for u in generators_U_plus] #write generators of U^+ in terms of powers of generators of U_K
        #we now use a free group to compute the quotient group
        F.<a,b> = FreeGroup()
        commutators = [F([1, 2, -1, -2])] #unit group is commutative
        roots_of_unity = [a^UK.zeta_order()] #first generator is a root of unity
        pos_units = [prod((u**e for u,e in zip([a,b],logs[0])), F(1) )] #write generators of U^+ into the free group
        relations = commutators + pos_units + roots_of_unity
        G = F/relations #quotient group
        elements = [el.Tietze() for el in G.list()] #output elements of quotient group in readable form
        UK.inject_variables(verbose=False)
        unit_list = [self.unit_maker(element) for element in elements]
        self._units_for_xi=unit_list
        return self._units_for_xi


    def good_generator(self,CM_type,princ_ideal):
        """
        Given a CMField K, a primitive CM type and PRINCIPAL ideal, returns b, a generator of the ideal that is totally imaginary and that has imaginary part negative under each embedding in the CM type
        """
        prec = self._prec

        delta = princ_ideal.gens_reduced()[0]

        try:
            units_for_xi = self._units_for_xi
        except:
            units_for_xi = self.units_for_xi()

        for u in units_for_xi:
            if all([CM_type[i](u*delta).imag()<0 for i in range(2)]) and (u*delta) + (u*delta).conjugate() == 0:
                return u*delta
        return None


    def princ_polarized(self,CM_type):
        """
        Given a CMField K and a CM-type Phi, returns a list of all pairs (ideal, xi) such that CC^3/Phi(ideal) is a principally polarized abelian variety with p.p. given by xi.
        """
        all_pairs = []
        conj = self.complex_conjugation()
        try:
            suitable_ideals = self._suitable_ideals
        except:
            suitable_ideals = self.suitable_ideals()
        for ideal in suitable_ideals:
            princ_ideal = self.different() * ideal * conj(ideal)
            b = self.good_generator(CM_type,princ_ideal)
            if not(b == None):
                xi = b**(-1)
                for u in self.units_epsilon():
                    all_pairs.append([ideal,u*xi])
        return all_pairs


    def xi_and_CMtype(self):
        allphis = self.all_CMtypes()
        pair=[]
        for phis in allphis:
            delta=K.good_generator(phis,principal)
            if not (delta == None):
                pair.append([delta, phis])
        return pair


def find_conj(embeddings,phi):
    """
    複素共役を探す
    """
    for embed in embeddings:
        if embed.im_gens()[0].conjugate() == phi.im_gens()[0]:
            return embed



def is_CMtype(phis):
    """
    CMタイプかどうか判定
    """
    if (phis[0].im_gens()[0].conjugate()==phis[1].im_gens()[0]):
        return False
    else:
        return True
