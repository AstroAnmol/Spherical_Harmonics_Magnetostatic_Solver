from spherical_harmonic_two_grain import spherical_harmonic_two_grain
import numpy as np
import pandas as pd

## Inputs
mu0 = 4*np.pi*1e-07
B0 = mu0
susc = 1.0 # Magnetic susceptibility
a = 1.0 # Grain radius, meters
sep=4    
alpha=0.0
L=10 # Number of multipoles used
debug_mag=1

Hr_L, Hth_L, Hphi_L, f_L=spherical_harmonic_two_grain(B0,susc, a, sep, alpha, L)

np.save('Hr_L_susc_1_L_10_sep_4', Hr_L)
np.save('Hth_L_susc_1_L_10_sep_4', Hth_L)
np.save('Hphi_L_susc_1_L_10_sep_4', Hphi_L)
np.save('f_L_susc_1_L_10_sep_4', f_L)

# fmag=np.zeros(12)

# sep=np.arange(2.0, 4.3, 0.2)
# check=range(12)

# for i in range(12):
#     print(i+1)
#     f=spherical_harmonic_two_grain(B0,susc, a, sep[i], alpha, L)
#     fmag[i]=f[2]
#     print(f[2])

# data_d={'Force':fmag, 'Separation':sep}
# data=pd.DataFrame(data=data_d)
# data.to_csv('fmag_susc_40_th_90_L_50.csv')

# fmag=np.zeros(5)

# a=np.arange(10.0, 60.0, 10.0)
# check=range(5)
# # print(a)
# # print(check)

# for i in range(5):
#     print(i+1)
#     f=spherical_harmonic_two_grain(B0,susc, a[i], sep, alpha, L)
#     fmag[i]=f[2]
#     print(f[2])

# data_d={'Force':fmag, 'particle_radius':a}
# # data_d={'Force':fmag, 'Separation':sep}
# data=pd.DataFrame(data=data_d)
# data.to_csv('fmag_susc_1_th_90_L_10_sep_2.csv')