import numpy as np
from scipy.special import comb
from scipy.special import lpmn
import matplotlib.pyplot as plt

## Inputs
mu0 = 4*np.pi*1e-07
B0 = mu0
susc = 30.0 # Magnetic susceptibility
a = 1.0 # Grain radius, meters
sep=4.0    
alpha=0.0
L=50 # Number of multipoles used
mu = (a+susc)*mu0
H0mag = B0/mu0 #Applied magnetic field, A/m
sep = sep*a

R1=np.array([[0], [0], [sep]])
R2=np.array([[0], [0], [0]])

alpha= np.deg2rad(alpha) #change angle into radians
H0= np.array([[H0mag*np.sin(alpha)], [0], [H0mag*np.cos(alpha)]]) # A/m


#Create a 3D spherical mesh
dang= np.pi/180.0
inc= np.arange(dang/2, np.pi + dang, dang)
az= np.arange(dang/2, 2*np.pi + dang, dang)
dr= a/100.0
r1= np.arange(a-dr,a+50*dr,dr)
theta, phi, R=np.meshgrid(inc,az,r1)
print(r1)

#define size parameters
size_R=np.shape(R) #size of all elements in 3D space
size_H=np.append(size_R,L) #size for 4D mag arrays

#Convert the spherical mesh to cartesian
# x= np.multiply(R,np.multiply(np.cos(phi),np.sin(theta)))
# y= np.multiply(R,np.multiply(np.sin(phi),np.sin(theta)))
# z= np.multiply(R,np.cos(theta))

l=L

x=R*np.cos(phi)*np.sin(theta)
y=R*np.sin(phi)*np.sin(theta)
z=R*np.cos(theta)

#Load mag fields
Hr_L=np.load('Hr_L.dat.npy')
Hth_L=np.load('Hth_L.dat.npy')
Hphi_L=np.load('Hphi_L.dat.npy')
f_L=np.load('f_L.dat.npy')

H_tot_mag=np.zeros(size_R)
for i in range(size_R[0]):
    for j in range(size_R[1]):
        for k in range(size_R[2]):
            ph=az[i]
            th=inc[j]
            #Transformation matrix
            pre= np.array([[np.sin(th)*np.cos(ph), np.cos(th)*np.cos(ph), -np.sin(ph)],
                    [np.sin(th)*np.sin(ph), np.cos(th)*np.sin(ph),  np.cos(ph)],
                    [np.cos(th), -np.sin(th),  0]])
            post=np.transpose(pre)
            H0_sph=np.matmul(post,H0)
            H_sph=np.array([[Hr_L[i,j,k,l-1]],[Hth_L[i,j,k,l-1]],[Hphi_L[i,j,k,l-1]]]) + H0_sph
            H_tot_mag[i,j,k]=np.linalg.norm(H_sph)

print(size_R)
plt.figure()
Ang=45
pc=plt.pcolor(np.squeeze(z[Ang,:,:])/a,np.squeeze(x[Ang,:,:])/a,np.squeeze(H_tot_mag[Ang,:,:]))
# plt.colormap('hot')
plt.xlim((-2, 2))
plt.ylim((-1.5, 1.5))
plt.title('|H|')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.axis('equal')
plt.grid(True)
plt.show()

# print(f_L)