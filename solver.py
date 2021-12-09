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
L=40 # Number of multipoles used
debug_mag=1

# f=spherical_harmonic_two_grain(B0,susc, a, sep, alpha, L)

# print(f)
fmag_0_1=np.zeros(12)

sep=np.arange(2.0, 4.4, 0.2)
susc=12
for i in range(13):
    print(i+1)
    f=spherical_harmonic_two_grain(B0,susc, a, sep[i], alpha, L)
    fmag_0_1[i]=f[2]
    print(f[2])
