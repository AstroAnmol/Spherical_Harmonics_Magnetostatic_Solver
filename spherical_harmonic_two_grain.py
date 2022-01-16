import numpy as np
from scipy.special import comb
from scipy.special import lpmn
# import matplotlib.pyplot as plt

@np.vectorize
def lpmn_arr(m, n, x):
    return lpmn(m, n, x)[0][-1, -1]
@np.vectorize
def d_lpmn_arr(m, n, x):
    return lpmn(m, n, x)[1][-1, -1]

def spherical_harmonic_two_grain(B0, susc, a, sep, alpha, L):
    mu0= 4*np.pi*1e-07
    mu = (a+susc)*mu0
    H0mag = B0/mu0 #Applied magnetic field, A/m
    sep = sep*a

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

    ## Computing Magnetic Field

    #Create a 3D spherical mesh
    dang= np.pi/180.0
    inc= np.arange(dang/2, np.pi + dang, dang)
    az= np.arange(dang/2, 2*np.pi + dang, dang)
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
        for m in range(2):
            Hrs=0
            Hths=0
            for s in np.arange(m,L+1):
                Psm=lpmn_arr(m,s,np.cos(theta))
                dPsm=d_lpmn_arr(m,s,np.cos(theta))*(-np.sin(theta))
                # R component
                additional=((-1)**(s+m))*s*np.multiply(np.power(R,s-1),Psm)/ (sep**(l+s+1))
                Hrs=Hrs + additional*comb(l+s,s+m, exact=True)
                # Theta component
                Hths=Hths + ((-1)**(s+m)) * comb(l+s,s+m, exact=True)*np.multiply(np.power(R,(s-1)),dPsm)/((sep**(l+s+1)))

            Plm=lpmn_arr(m,l,np.cos(theta))
            dPlm=d_lpmn_arr(m,l,np.cos(theta))*(-np.sin(theta))
            if m==0:
                # R Component
                Hr=Hr + np.multiply(((l+1)*Beta1_0[l-1]*np.divide(Plm,np.power(R,(l+2))) -  Beta2_0[l-1]*Hrs),np.cos(m*phi))
                # Theta component
                Hth=Hth + np.multiply((Beta1_0[l-1]*np.divide(dPlm,np.power(R,(l+2))) + Beta2_0[l-1]*Hths),np.cos(m*phi))
            elif m==1:
                # R Component
                Hr=Hr + np.multiply(((l+1)*Beta1_1[l-1]*np.divide(Plm,np.power(R,(l+2))) -  Beta2_1[l-1]*Hrs),np.cos(m*phi))
                # Theta component
                Hth=Hth + np.multiply((Beta1_1[l-1]*np.divide(dPlm,np.power(R,(l+2))) + Beta2_1[l-1]*Hths),np.cos(m*phi))
        Hr_L[:,:,:,l-1]=Hr
        Hth_L[:,:,:,l-1]=-Hth

    # Phi compenent
    for l in np.arange(1,L+1):
        Hphis=0
        for s in np.arange(1,L+1):
            Ps1=lpmn_arr(1,s,np.cos(theta))
            Hphis=Hphis + (-1)**(s+1) *comb(l+s,s+1, exact=True)*np.divide(np.multiply(np.power(R,s-1),Ps1),np.sin(theta))/(sep**(l+s+1))
        Pl1=lpmn_arr(1,l,np.cos(theta))
        # Hphi=Hphi + np.multiply((Beta1_1[l-1]*np.divide(Pl1,np.multiply(np.sin(theta),np.power(R,l+2))) + Beta2_1[l-1]*Hphis),np.sin(phi))
        Hphi_L[:,:,:,l-1]=Hphi

    Hth=-Hth

    # Formulating the Maxwell Stress Tensor in Spherical Coordinates
    f_L=np.zeros((3,1,L))
    for l in np.arange(1,L+1):
        f=0
        for i in range(size_R[0]):
            if i==0 or i==(size_R[0]-1):
                p=1
            elif np.mod(i,2)==0:
                p=2
            else:
                p=4
            for j in range(size_R[1]):
                if j==0 or j==size_R[1]-1:
                    q=1
                elif np.mod(j,2)==0:
                    q=2
                else:
                    q=4
                ph=az[i]
                th=inc[j]
                # Transformation matrix
                pre= np.array([[np.sin(th)*np.cos(ph), np.cos(th)*np.cos(ph), -np.sin(ph)],
                                [np.sin(th)*np.sin(ph), np.cos(th)*np.sin(ph),  np.cos(ph)],
                                [np.cos(th), -np.sin(th),  0]])
                post=np.transpose(pre)
                H0_sph=np.matmul(post,H0)
                H_sph=np.array([[Hr_L[i,j,1,l-1]],[Hth_L[i,j,1,l-1]],[Hphi_L[i,j,1,l-1]]]) + H0_sph
                H_cart=np.matmul(pre,H_sph)
                h=np.linalg.norm(H_cart)
                T_cart=mu0*(np.matmul(H_cart,np.transpose(H_cart)) - 0.5*(h**2)*np.eye(3))
                rn_hat=np.array([[np.sin(th)*np.cos(ph)],[np.sin(th)*np.sin(ph)],[np.cos(th)]])
                f=f+ a*a*np.sin(th)*p*q*np.matmul(T_cart,rn_hat)
        f=f*dang*dang/9.0
        f_L[:,:,l-1]=f  
    return(f)
# #plot z component of force vs L
# if debug_f_L==1:
#     plt.figure()
#     plt.plot(np.arange(1,L+1),f_L[2,0,:])
#     plt.title("Z component of force vs L (multipoles)")
#     plt.xlabel('L (multipoles)')
#     plt.ylabel('f_z (N)')
#     plt.grid(True)
#     plt.show()