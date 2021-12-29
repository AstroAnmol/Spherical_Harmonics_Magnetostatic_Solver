from spherical_harmonic_two_grain import spherical_harmonic_two_grain
import numpy as np
import pandas as pd

## Inputs
mu0 = 4*np.pi*1e-07
B0 = mu0
susc = 1.0 # Magnetic susceptibility
a = 1.0 # Grain radius, meters
sep=2.0    
alpha=0.0
L=80 # Number of multipoles used
debug_mag=1

# f=spherical_harmonic_two_grain(B0,susc, a, sep, alpha, L)

# print(f)
fmag=np.zeros(12)

sep=np.arange(2.0, 4.3, 0.2)
susc=75
check=range(12)

for i in range(12):
    print(i+1)
    f=spherical_harmonic_two_grain(B0,susc, a, sep[i], alpha, L)
    fmag[i]=f[2]
    print(f[2])

data_d={'Force':fmag, 'Separation':sep}
data=pd.DataFrame(data=data_d)
data.to_csv('fmag_susc_75_th_0_L_80.csv')