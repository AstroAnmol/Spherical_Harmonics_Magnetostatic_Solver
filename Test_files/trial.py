import numpy as np
from scipy.special import comb
from scipy.special import lpmn
import time


@np.vectorize
def lpmn_arr(m, n, x):
    return lpmn(m, n, x)[0][-1, -1]
@np.vectorize
def d_lpmn_arr(m, n, x):
    return lpmn(m, n, x)[1][-1, -1]


start_time = time.time()


## Inputs
mu0 = 4*np.pi*1e-07
B0 = mu0
susc = 1.0 # Magnetic susceptibility
a = 1.0 # Grain radius, meters
sep=2    
alpha=0.0
L=10 # Number of multipoles used
debug_mag=0
mu = (1+susc)*mu0
H0mag = B0/mu0 #Applied magnetic field, A/m
sep = sep*a
debug_mag=0

R1=np.array([[0], [0], [sep]])
R2=np.array([[0], [0], [0]])

alpha= np.deg2rad(alpha) #change angle into radians
H0= np.array([[H0mag*np.sin(alpha)], [0], [H0mag*np.cos(alpha)]]) # A/m

# parallel and perpendicular components of the applied magnetic field
H_perp=H0[0]
H_prll=H0[2]

# Code to solve for coefficents (using paper)
for m in range(2):
    # Creating the L X L matrices
    X= np.zeros((L,L))
    Delta_m= np.zeros((L,L))
    Gamma_m= np.zeros((L,L))
    for i in range(L):
        for j in range(L):
            # X Matrix
            if i==j:
                X[i,j]=(i+1)*(mu/mu0) + (i+1) + 1
            # Delta and Gamma matrix
            Delta_m[i,j]= ((-1)**((i+1)+m))*((i+1)*(mu/mu0)-(i+1))*comb(i+1+j+1, j+1-m, exact=True)*(a**(2*(i+1)+1))/(sep**(i+1+j+1+1))
            Gamma_m[i,j]= ((-1)**(i+1+j+1))*Delta_m[i,j]
            # 2L X 2L matrix
    Am=np.zeros((2*L,2*L))
    Am[0:L,0:L]=X
    Am[L:,0:L]=Gamma_m
    Am[0:L,L:]=Delta_m
    Am[L:,L:]=X

    # q vector
    qm= np.zeros((L,1))
    if m==0:
        qm[0,0]=-H_prll*(a**3)*(1-mu/mu0)
    elif m==1:
        qm[0,0]=H_perp*(a**3)*(1-mu/mu0)
        
    # 2L Q vector
    Qm=np.append(qm, qm)

    #solve linear equations
    Beta_m=np.linalg.solve(Am,Qm)
    Beta1_m=Beta_m[0:L]
    Beta2_m=Beta_m[L:]
    if m==0:
        Beta1_0=Beta1_m
        Beta2_0=Beta2_m
    elif m==1:
        Beta1_1=Beta1_m
        Beta2_1=Beta2_m

print("Linear system solved")

def mag_field(R, theta, phi):
    Hr=0
    Hth=0
    Hphi=0
    for l in np.arange(1,L+1):
        for m in range(2):
            Hrs=0
            Hths=0
            Hphis=0
            for s in np.arange(m,L+1):
                Psm=lpmn_arr(m,s,np.cos(theta))
                dPsm=d_lpmn_arr(m,s,np.cos(theta))*(-np.sin(theta))
                ls_choose_sm=comb(l+s,s+m, exact=True)
                R_power_s1=np.power(R,(s-1))
                sep_power=sep**(l+s+1)
                minus_one_power=(-1)**(s+m)
                R_power_times_Psm=np.multiply(R_power_s1,Psm)

                # R component
                additional=(minus_one_power)*s*R_power_times_Psm/(sep_power)
                Hrs=Hrs + additional*ls_choose_sm
                # Theta component
                Hths=Hths + (minus_one_power)*ls_choose_sm*np.multiply(R_power_s1,dPsm)/(sep_power)
                # Phi component
                if m==1:
                    Hphis=Hphis + (minus_one_power)*ls_choose_sm*np.divide(R_power_times_Psm,np.sin(theta))/(sep_power)

            Plm=lpmn_arr(m,l,np.cos(theta))
            dPlm=d_lpmn_arr(m,l,np.cos(theta))*(-np.sin(theta))
            R_power_l2=np.power(R,(l+2))


            if m==0:
                # R Component
                Hr=Hr + np.multiply(((l+1)*Beta1_0[l-1]*np.divide(Plm,R_power_l2) -  Beta2_0[l-1]*Hrs),np.cos(m*phi))
                # Theta component
                Hth=Hth + np.multiply((Beta1_0[l-1]*np.divide(dPlm,R_power_l2) + Beta2_0[l-1]*Hths),np.cos(m*phi))
            elif m==1:
                # R Component
                Hr=Hr + np.multiply(((l+1)*Beta1_1[l-1]*np.divide(Plm,R_power_l2) -  Beta2_1[l-1]*Hrs),np.cos(m*phi))
                # Theta component
                Hth=Hth + np.multiply((Beta1_1[l-1]*np.divide(dPlm,R_power_l2) + Beta2_1[l-1]*Hths),np.cos(m*phi))
                # Phi component
                Hphi=Hphi + np.multiply((Beta1_1[l-1]*np.divide(Plm,np.multiply(np.sin(theta),R_power_l2)) + Beta2_1[l-1]*Hphis),np.sin(phi))


    Hth=-Hth

    return np.array([[Hr], [Hth], [Hphi]])


def integrand(th, ph):
    # Transformation matrix
    pre= np.array([[np.sin(th)*np.cos(ph), np.cos(th)*np.cos(ph), -np.sin(ph)],
                    [np.sin(th)*np.sin(ph), np.cos(th)*np.sin(ph),  np.cos(ph)],
                    [np.cos(th), -np.sin(th),  0]])
    post=np.transpose(pre)
    H0_sph=np.matmul(post,H0)
    H_sph= mag_field(a, th, ph) + H0_sph
    H_cart=np.matmul(pre,H_sph)
    h=np.linalg.norm(H_cart)
    T_cart=mu0*(np.matmul(H_cart,np.transpose(H_cart)) - 0.5*(h**2)*np.eye(3))
    rn_hat=np.array([[np.sin(th)*np.cos(ph)],[np.sin(th)*np.sin(ph)],[np.cos(th)]])
    return np.sin(th)*np.matmul(T_cart,rn_hat)

# def f_x(th, ph):
#     return integrand(th, ph)[0]

# force=dblquad(f_x, 0, np.pi, lambda ph: 0, lambda ph: 2*np.pi)


#Create a 3D spherical mesh
dang= np.pi/180
inc= np.arange(dang/2, np.pi + dang, dang)
az= np.arange(dang/2, 2*np.pi + dang, dang)

size_az=np.shape(az)[0]
size_inc=np.shape(inc)[0] #size of all elements in 3D space

# Formulating the Maxwell Stress Tensor in Spherical Coordinates
f=0
for i in range(size_az):
    if i==0 or i==(size_az-1):
        p=1
    elif np.mod(i,2)==0:
        p=2
    else:
        p=4
    for j in range(size_inc):
        if j==0 or j==size_inc-1:
            q=1
        elif np.mod(j,2)==0:
            q=2
        else:
            q=4
        ph=az[i]
        th=inc[j]
        # Transformation matrix
        f=f+ a*p*q*integrand(th, ph)
        # print(i,j)
f=f*dang*dang/9.0

print("--- %s seconds ---" % (time.time() - start_time))

print(f)