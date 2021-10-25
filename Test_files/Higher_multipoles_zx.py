import numpy as np
from scipy.special import comb
from scipy.special import lpmn
# def spherical_harmonic_two_grain(B0, susc, a, sep, alpha, L, debug_mag):

mu0 = 4*np.pi*1e-07
B0 = mu0 
susc = 1  
a = 1      
sep=2   
alpha=0 
L=10   

mu0 = 4*np.pi*1e-07
mu = (1+susc)*mu0
H0mag = B0/(mu0)  # Applied magnetic field, A/m
sep=sep*a

R1=np.array([0, 0, sep])
R2=np.array([0, 0, 0])

alpha = np.deg2rad(alpha) #change angle into radians
H0 = np.array([H0mag*np.sin(alpha), 0, H0mag*np.cos(alpha)]) # A/m

# parallel and perpendicular components of the applied magnetic field
H_perp=H0[0]
H_prll=H0[2]

# Code to solve for coefficients (using paper)
for m in range(1):
    # Creating the L X L matrices
    X=np.zeros((L,L))
    Delta_m=np.zeros((L,L))
    Gamma_m=np.zeros((L,L))
    for i in range(L):
        for j in range(L):
            # X Matrix
            if i==j:
                X[i,j]=(i+1)*(mu/mu0) + (i+1) + 1
        # Delta and Gamma matrix
        Delta_m[i,j]= ((-1)**((i+1)+m))*((i+1)*(mu/mu0)-(i+1))*comb(i+1+j+1, j+1-m)*(a^(2*(i+1)+1))/(sep^(i+1+j+1+1))
        Gamma_m[i,j]= ((-1)**(i+1+j+1))
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

## Computing Magnetic Field

#Create a 3D spherical mesh
dang= np.pi/18.0
inc= np.arange(dang/2, np.pi + dang/2, dang)
az= np.arange(dang/2, 2*np.pi + dang/2, dang)
dr= a/100.0
r1= np.arange(a-dr,a+dr,dr)
theta, phi, R=np.meshgrid(inc,az,r1)

#Convert the spherical mesh to cartesian
x= np.multiply(R,np.multiply(np.cos(phi),np.sin(theta)))
y= np.multiply(R,np.multiply(np.sin(phi),np.sin(theta)))
z= np.multiply(R,np.cos(theta))

#define size parameters
size_R=np.shape(R) #size of all elements in 3D space
size_H=np.append(size_R,L) #size for 4D mag arrays

Hr_L=np.zeros(size_H)
Hth_L=np.zeros(size_H)
Hphi_L=np.zeros(size_H)

Hr=0
Hth=0
Hphi=0

for l in np.arange(1,L+1):
    for m in range(1):
        Hrs=0
        Hths=0
        # for s in np.arange(m,L+1):

print(lpmn(10,10,np.array([0.1,0.2,0.3]))[0])

