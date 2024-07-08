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
drainage = 'undrained'
# number of cycles
Ncycles = 10
# for stress controlled defined CSR > 0
CSR = 0.15
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
t.CyDSS(
    sig3,
    'constant volume', 
    k0, 
    Ncycles, 
    Epslim, 
    CSR,
    figname=r'.\output\CyDSS_res.png',
    filename=r'.\output\CyDSS_res.csv')
