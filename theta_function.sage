"""

 This file has the functions to compute values of theta functions

"""
import warnings
warnings.simplefilter("ignore", UserWarning) #In order to remove warnings when computing eigenvalues


def theta_function(period_matrix, vec1, vec2, prec=664):
    a= period_matrix[0][0]
    b= period_matrix[0][1]
    c= period_matrix[1][1]

    eig = min((period_matrix.apply_map(lambda aux: aux.imag_part())).eigenvalues())

    dig = floor(RR(prec* log(2,10)))

    bound = ceil(RR(1./2 + sqrt(1./4 + dig*ln(10)/(pi.n(prec)*eig))))+10

    S =  'A = %s;'%a
    S += 'B = %s;'%b
    S += 'C = %s;'%c
    S += 'bound = %s;' %bound
    S += 'dig = %s;' %dig
    S +='default(realprecision,dig); \
    thetaval(del,eps,a,b,c,bd)= \
    {s=0; t=0; \
    cutoff = -dig*log(10.0)-2*log(bd); \
    for(i=-bd,bd, \
        for(j=-bd,bd, \
            t = Pi*I*((i+1/2*del[1])^2*a + 2*(i+1/2*del[1])*(j+1/2*del[2])*b +(j+1/2*del[2])^2*c) \
                +2*Pi*I*(1/2* (i + 1/2*del[1])*eps[1]+1/2*(j+1/2*del[2])*eps[2]);\
    if(real(t)>cutoff,s = s + exp(t))));  return(s);}'

    gp(S)
    V =  'v1 = %s;'%vec1[0]
    V += 'v2 = %s;'%vec1[1]
    V += 'v3 = %s;'%vec2[0]
    V += 'v4 = %s;'%vec2[1]
    V+= 'Vec1=[v1,v2];'
    V+= 'Vec2=[v3,v4];'
    V+= 'thetan=thetaval(Vec1,Vec2,A,B,C,bound);'
    gp(V)
    theta=gp.eval('thetan')
    Cec=ComplexField(prec)
    return Cec(theta)

def theta_constants(period_matrix,prec=664):
    all_evens = [[[0,0],[0,0]],[[0,0],[1,0]],[[0,0],[0,1]],[[0,0],[1,1]],[[1,0],[0,0]],[[1,0],[0,1]], \
                [[0,1],[0,0]],[[0,1],[1,0]],[[1,1],[0,0]],[[1,1],[1,1]]]
    all_values = [theta_function(period_matrix,even[0],even[1],prec) for even in all_evens]
    return all_values

def h10(period_matrix,prec=664):
    t1,t2,t3,t4,t5,t6,t7,t8,t9,t10=theta_constants(period_matrix,prec)
    S = 'T1 = %s;'%t1
    S +='T2 = %s;'%t2
    S +='T3 = %s;'%t3
    S +='T4 = %s;'%t4
    S +='T5 = %s;'%t5
    S +='T6 = %s;'%t6
    S +='T7 = %s;'%t7
    S +='T8 = %s;'%t8
    S +='T9 = %s;'%t9
    S +='T10 = %s;'%t10
    S +='s = T1^2*T2^2*T3^2*T4^2*T5^2*T6^2*T7^2*T8^2*T9^2*T10^2;'
    gp(S)
    h_10=gp.eval('s')
    Cec=ComplexField(prec)
    return Cec(h_10)

def h16(period_matrix,prec=664):
    t1,t2,t3,t4,t5,t6,t7,t8,t9,t10=theta_constants(period_matrix,prec)
    S = 't1 = %s;'%t1
    S +='t2 = %s;'%t2
    S +='t3 = %s;'%t3
    S +='t4 = %s;'%t4
    S +='t5 = %s;'%t5
    S +='t6 = %s;'%t6
    S +='t7 = %s;'%t7
    S +='t8 = %s;'%t8
    S +='t9 = %s;'%t9
    S +='t10 = %s;'%t10
    S +='s = t8^4*t1^4*t5^4*t2^4*t9^4*t6^4*t8^4*t10^4+t5^4*t1^4*t5^4*t2^4*t9^4*t6^4*t8^4*t3^4+t10^4*t1^4*t2^4*t9^4*t6^4*t8^4*t10^4*t3^4+\
           t3^4*t1^4*t5^4*t2^4*t9^4*t6^4*t10^4*t3^4+t1^4*t1^4*t5^4*t9^4*t6^4*t8^4*t10^4*t7^4+t2^4*t5^4*t2^4*t9^4*t6^4*t8^4*t10^4*t7^4+\
           t1^4*t1^4*t5^4*t2^4*t6^4*t8^4*t3^4*t7^4+t9^4*t5^4*t2^4*t9^4*t6^4*t8^4*t3^4*t7^4+t9^4*t1^4*t5^4*t2^4*t9^4*t6^4*t8^4*t10^4*t7^4+\
           t6^4*t1^4*t5^4*t2^4*t6^4*t10^4*t3^4*t7^4+t5^4*t1^4*t5^4*t9^4*t8^4*t10^4*t3^4*t7^4+t2^4*t1^4*t2^4*t9^4*t8^4*t10^4*t3^4*t7^4+\
           t6^4*t1^4*t9^4*t6^4*t8^4*t10^4*t3^4*t7^4+t8^4*t1^4*t5^4*t2^4*t8^4*t10^4*t3^4*t7^4+t10^4*t5^4*t2^4*t6^4*t8^4*t10^4*t3^4*t7^4+\
           t3^4*t5^4*t9^4*t6^4*t8^4*t10^4*t3^4*t7^4 +t7^4*t1^4*t5^4*t2^4*t9^4*t6^4*t10^4*t7^4+t7^4*t1^4*t2^4*t9^4*t6^4*t8^4*t3^4*t7^4+\
           t9^4*t1^4*t5^4*t2^4*t9^4*t8^4*t10^4*t4^4+t6^4*t1^4*t5^4*t2^4*t6^4*t8^4*t10^4*t4^4+t2^4*t1^4*t5^4*t2^4*t9^4*t8^4*t3^4*t4^4+\
           t6^4*t1^4*t5^4*t9^4*t6^4*t8^4*t3^4*t4^4+t1^4*t1^4*t5^4*t9^4*t6^4*t10^4*t3^4*t4^4+t2^4*t5^4*t2^4*t9^4*t6^4*t10^4*t3^4*t4^4+\
           t1^4*t1^4*t2^4*t6^4*t8^4*t10^4*t3^4*t4^4+t5^4*t5^4*t2^4*t6^4*t8^4*t10^4*t3^4*t4^4+t9^4*t2^4*t9^4*t6^4*t8^4*t10^4*t3^4*t4^4+\
           t8^4*t5^4*t9^4*t6^4*t8^4*t10^4*t3^4*t4^4+t10^4*t1^4*t5^4*t9^4*t8^4*t10^4*t3^4*t4^4+t3^4*t1^4*t5^4*t2^4*t8^4*t10^4*t3^4*t4^4+\
           t5^4*t1^4*t5^4*t2^4*t9^4*t6^4*t7^4*t4^4+t2^4*t1^4*t5^4*t2^4*t6^4*t8^4*t7^4*t4^4+t9^4*t1^4*t5^4*t9^4*t6^4*t8^4*t7^4*t4^4+\
           t8^4*t1^4*t2^4*t9^4*t6^4*t8^4*t7^4*t4^4+t1^4*t1^4*t2^4*t9^4*t8^4*t10^4*t7^4*t4^4+t5^4*t5^4*t2^4*t9^4*t8^4*t10^4*t7^4*t4^4+\
           t6^4*t2^4*t9^4*t6^4*t8^4*t10^4*t7^4*t4^4+t10^4*t1^4*t2^4*t9^4*t6^4*t10^4*t7^4*t4^4+t10^4*t1^4*t5^4*t6^4*t8^4*t10^4*t7^4*t4^4+\
           t1^4*t1^4*t5^4*t2^4*t9^4*t3^4*t7^4*t4^4+t6^4*t5^4*t2^4*t9^4*t6^4*t3^4*t7^4*t4^4+t8^4*t5^4*t2^4*t9^4*t8^4*t3^4*t7^4*t4^4+\
           t5^4*t1^4*t5^4*t6^4*t10^4*t3^4*t7^4*t4^4+t2^4*t1^4*t2^4*t6^4*t10^4*t3^4*t7^4*t4^4+t9^4*t1^4*t9^4*t6^4*t10^4*t3^4*t7^4*t4^4+\
           t8^4*t1^4*t6^4*t8^4*t10^4*t3^4*t7^4*t4^4+t10^4*t5^4*t2^4*t9^4*t10^4*t3^4*t7^4*t4^4+t3^4*t1^4*t2^4*t9^4*t6^4*t3^4*t7^4*t4^4+\
           t3^4*t1^4*t5^4*t6^4*t8^4*t3^4*t7^4*t4^4+t3^4*t2^4*t9^4*t8^4*t10^4*t3^4*t7^4*t4^4+t7^4*t1^4*t5^4*t2^4*t8^4*t10^4*t7^4*t4^4+\
           t7^4*t1^4*t5^4*t9^4*t8^4*t3^4*t7^4*t4^4+t7^4*t5^4*t9^4*t6^4*t10^4*t3^4*t7^4*t4^4+t7^4*t2^4*t6^4*t8^4*t10^4*t3^4*t7^4*t4^4+\
           t4^4*t1^4*t5^4*t2^4*t9^4*t6^4*t10^4*t4^4+t4^4*t1^4*t2^4*t9^4*t6^4*t8^4*t3^4*t4^4+t4^4*t5^4*t9^4*t6^4*t8^4*t10^4*t7^4*t4^4+\
           t4^4*t5^4*t2^4*t6^4*t8^4*t3^4*t7^4*t4^4+t4^4*t1^4*t5^4*t2^4*t10^4*t3^4*t7^4*t4^4+t4^4*t1^4*t9^4*t8^4*t10^4*t3^4*t7^4*t4^4;'
    gp(S)
    h_16=gp.eval('s')
    Cec=ComplexField(prec)
    return Cec(h_16)

def Igusa_Clebsch(period_matrix,prec=664):
    t1,t2,t3,t4,t5,t6,t7,t8,t9,t10=theta_constants(period_matrix,prec)
    h_4=0
    for i in range(10):
        x = theta_constants(period_matrix,prec)[i]^8
        h_4 = h_4 + x

    h_10 = s=t1^2*t2^2*t3^2*t4^2*t5^2*t6^2*t7^2*t8^2*t9^2*t10^2
    h_12 = (t1*t5*t2*t9*t6*t10)^4+(t1*t2*t9*t6*t8*t3)^4+(t5*t9*t6*t8*t10*t7)^4+(t5*t2*t6*t8*t3*t7)^4+\
           (t1*t5*t2*t10*t3*t7)^4+ (t1*t9*t8*t10*t3*t7)^4+(t1*t5*t2*t8*t10*t4)^4+(t1*t5*t9*t8*t3*t4)^4 +\
           (t5*t9*t6*t10*t3*t4)^4+(t2*t6*t8*t10*t3*t4)^4+(t1*t2*t9*t6*t7*t4)^4+(t1*t5*t6*t8*t7*t4)^4+\
           (t2*t9*t8*t10*t7*t4)^4+(t5*t2*t9*t3*t7*t4)^4+(t1*t6*t10*t3*t7*t4)^4

    h_16 = t8^4*(t1*t5*t2*t9*t6*t8*t10)^4+t5^4*(t1*t5*t2*t9*t6*t8*t3)^4+t10^4*(t1*t2*t9*t6*t8*t10*t3)^4+\
           t3^4*(t1*t5*t2*t9*t6*t10*t3)^4+t1^4*(t1*t5*t9*t6*t8*t10*t7)^4+t2^4*(t5*t2*t9*t6*t8*t10*t7)^4+\
           t1^4*(t1*t5*t2*t6*t8*t3*t7)^4+t9^4*(t5*t2*t9*t6*t8*t3*t7)^4+t9^4*(t1*t5*t2*t9*t10*t3*t7)^4+\
           t6^4*(t1*t5*t2*t6*t10*t3*t7)^4+t5^4*(t1*t5*t9*t8*t10*t3*t7)^4+t2^4*(t1*t2*t9*t8*t10*t3*t7)^4+\
           t6^4*(t1*t9*t6*t8*t10*t3*t7)^4+t8^4*(t1*t5*t2*t8*t10*t3*t7)^4+t10^4*(t5*t2*t6*t8*t10*t3*t7)^4+\
           t3^4*(t5*t9*t6*t8*t10*t3*t7)^4 +t7^4*(t1*t5*t2*t9*t6*t10*t7)^4+t7^4*(t1*t2*t9*t6*t8*t3*t7)^4+\
           t9^4*(t1*t5*t2*t9*t8*t10*t4)^4+t6^4*(t1*t5*t2*t6*t8*t10*t4)^4+t2^4*(t1*t5*t2*t9*t8*t3*t4)^4+\
           t6^4*(t1*t5*t9*t6*t8*t3*t4)^4+t1^4*(t1*t5*t9*t6*t10*t3*t4)^4+t2^4*(t5*t2*t9*t6*t10*t3*t4)^4+\
           t1^4*(t1*t2*t6*t8*t10*t3*t4)^4+t5^4*(t5*t2*t6*t8*t10*t3*t4)^4+t9^4*(t2*t9*t6*t8*t10*t3*t4)^4+\
           t8^4*(t5*t9*t6*t8*t10*t3*t4)^4+t10^4*(t1*t5*t9*t8*t10*t3*t4)^4+t3^4*(t1*t5*t2*t8*t10*t3*t4)^4+\
           t5^4*(t1*t5*t2*t9*t6*t7*t4)^4+t2^4*(t1*t5*t2*t6*t8*t7*t4)^4+t9^4*(t1*t5*t9*t6*t8*t7*t4)^4+\
           t8^4*(t1*t2*t9*t6*t8*t7*t4)^4+t1^4*(t1*t2*t9*t8*t10*t7*t4)^4+t5^4*(t5*t2*t9*t8*t10*t7*t4)^4+\
           t6^4*(t2*t9*t6*t8*t10*t7*t4)^4+t10^4*(t1*t2*t9*t6*t10*t7*t4)^4+t10^4*(t1*t5*t6*t8*t10*t7*t4)^4+\
           t1^4*(t1*t5*t2*t9*t3*t7*t4)^4+t6^4*(t5*t2*t9*t6*t3*t7*t4)^4+t8^4*(t5*t2*t9*t8*t3*t7*t4)^4+\
           t5^4*(t1*t5*t6*t10*t3*t7*t4)^4+t2^4*(t1*t2*t6*t10*t3*t7*t4)^4+t9^4*(t1*t9*t6*t10*t3*t7*t4)^4+\
           t8^4*(t1*t6*t8*t10*t3*t7*t4)^4+t10^4*(t5*t2*t9*t10*t3*t7*t4)^4+t3^4*(t1*t2*t9*t6*t3*t7*t4)^4+\
           t3^4*(t1*t5*t6*t8*t3*t7*t4)^4+t3^4*(t2*t9*t8*t10*t3*t7*t4)^4+t7^4*(t1*t5*t2*t8*t10*t7*t4)^4+\
           t7^4*(t1*t5*t9*t8*t3*t7*t4)^4+t7^4*(t5*t9*t6*t10*t3*t7*t4)^4+t7^4*(t2*t6*t8*t10*t3*t7*t4)^4+\
           t4^4*(t1*t5*t2*t9*t6*t10*t4)^4+t4^4*(t1*t2*t9*t6*t8*t3*t4)^4+t4^4*(t5*t9*t6*t8*t10*t7*t4)^4+\
           t4^4*(t5*t2*t6*t8*t3*t7*t4)^4+t4^4*(t1*t5*t2*t10*t3*t7*t4)^4+t4^4*(t1*t9*t8*t10*t3*t7*t4)^4


    I_2=h_12/h_10
    I_4=h_4
    I_6=h_16/h_10
    I_10=h_10

    return I_2,I_4,I_6,I_10


def Igusa_invariants(period_matrix,prec=664):
    A,B,C,D = Igusa_Clebsch(period_matrix,prec)

    j_1 = A^5/D
    j_2 = B*A^3/D
    j_3= C*A^2/D

    return j_1, j_2, j_3
