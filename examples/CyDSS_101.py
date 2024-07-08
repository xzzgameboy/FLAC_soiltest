'''
Copyright SRK Consulting (Canada) Inc. 2022 - All Rights Reserved
Licensed under CC BY-SA https://creativecommons.org/licenses/by-sa/4.0/
Author: info@srk.com

================================================================================
Cyclic Simple Direct Simple Shear Test

    - Drained, undrained, constant volume
    - Isotropic or anisotropic
    
================================================================================
'''

from src import SoilTest as st

t = st.SoilTest()

# ======================================================================== INPUT
# effective confining stress
sig3 = 100.
# effective stress ratio k0 = sig3/sig1
k0 = 1.0
# drainage conditions: 'drained', 'undrained', 'constant volume'
drainage = 'constant volume'
# number of cycles
Ncycles = 3
# for stress controlled defined CSR > 0
CSR = 0.10
# for strain controlled defined Epslim > 0
Epslim = 0.
# constitutive model and parameters
t.model = {
    'model': 'p2psand',
    'stress-xx-initial'  : -sig3,# Gamma when C3 = 0
    'stress-yy-initial'  : -sig3,# lambda_e when C3 = 0
    'stress-zz-initial'  : -sig3/k0,
    'stress-xy-initial': 0.,
    'stress-yz-initial': 0.,
    'stress-xz-initial': 0.,    
    'relative-density-initial'  : 0.35,# initial relative density
    }
    
# ===================================================================== run test
#Cyclic DSS test
t.CyDSS(sig3,drainage,k0,Ncycles=Ncycles,Epslim=0.0,CSR=CSR)
# Monotonic axial strain level
Epslim = 4.
t.calc_DSScv(epsyz_max=Epslim)
# save results
t.plot_ps (figname = r'.\output\DSS_res.png')
t.save_res (filename = r'.\output\DSS_res.csv')
