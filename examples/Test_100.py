'''
Copyright SRK Consulting (Canada) Inc. 2022 - All Rights Reserved
Licensed under CC BY-SA https://creativecommons.org/licenses/by-sa/4.0/
Author: info@srk.com

================================================================================
Simple Soil Test

    - Triaxial Compression - TXX
    - Direct Simple Shear - DSS
    - Oedometer - OED
    - Cyclic Triaxial - CyTXX
    - Cyclic DSS - CyDSS
    
================================================================================
'''

from src import SoilTest as st
t = st.SoilTest()

# ============================================ constitutive model and parameters
t.model = {
    'model': 'plastic-hardening',
    'over-consolidation-ratio': 1.0,
    'stiffness-ur-reference': 90.e3, 
    'stiffness-50-reference': 30.e3,
    'stiffness-oedometer-reference': 30.e3, 
    'pressure-reference' : 100.,
    'exponent' : 0.5,
    'friction' : 30.
    } 


# ============================================================ run triaxial test
# effective confining stress
sig3 = 70.
# effective stress ratio k0 = sig3/sig1
k0 = 1.0
# drainage conditions: 'drained', 'undrained', 'constant volume'
drainage = 'constant volume'
# maximum axial strain
epsa_max = 2.5
# Initial stress state parameters
t.model['stress-1-effective'] = -sig3/k0
t.model['stress-2-effective'] = -sig3
t.model['stress-3-effective'] = -sig3
# run test
t.TXX(
    sig3=sig3, 
    k0=k0, 
    drainage=drainage, 
    epsa_max = epsa_max,
    figname = r'.\output\TXX_res.png', 
    filename=r'.\output\TXX_res.csv')
    
    
'''
# ========================================================================== DSS
# effective confining stress
sig3 = 100.
# effective stress ratio k0 = sig3/sig1
k0 = 1.0
# drainage conditions: 'drained', 'undrained', 'constant volume'
drainage = 'constant volume'
# maximum shear strain
epsyz_max = 2.5
# Initial stress state parameters
t.model['stress-1-effective']= -sig3/k0
t.model['stress-2-effective']= -sig3
t.model['stress-3-effective']= -sig3

#run test
t.DSS(
    sig3, 
    k0, 
    drainage, 
    epsyz_max,
    figname = r'.\output\DSS_res.png', 
    filename= r'.\output\DSS.csv')
'''
''' 
# =========================================================== run oedometer test
# Initial stress state parameters
t.model['stress-1-effective'] = 0.
t.model['stress-2-effective'] = 0.
t.model['stress-3-effective'] = 0.
# load increments
load = [10, 100, 200, 400, 800, 400, 200, 400, 800, 1600, 3200, 1600, 800]

# run test
t.OED(
    load,
    figname = r'.\output\OED_res.png', 
    filename=r'.\output\OED_res.csv')
'''   
'''
# ============================================================== Cyclic Triaxial
# effective confining stress
sig3 = 70.
# effective stress ratio k0 = sig3/sig1
k0 = 1.0
# drainage conditions: 'drained', 'undrained', 'constant volume'
drainage = 'constant volume'
# number of cycles
Ncycles = 5
# for stress controlled defined CSR > 0
CSR = 0.1
# for strain controlled defined Epslim > 0
Epslim = 0.
# HSS parameters
t.model['stress-1-effective']= -sig3/k0
t.model['stress-2-effective']= -sig3
t.model['stress-3-effective']= -sig3

t.CyTXX(
    sig3, 
    drainage, 
    Ncycles=Ncycles, 
    Epslim=Epslim, 
    CSR=CSR,
    figname = r'.\output\CyTXX_res.png', 
    filename= r'.\output\CyTXX_res.csv')
'''
''' 
# =================================================================== Cyclic DSS
# effective confining stress
sig3 = 100.
# effective stress ratio k0 = sig3/sig1
k0 = 1.0
# drainage conditions: 'drained', 'undrained', 'constant volume'
drainage = 'undrained'
# number of cycles
Ncycles = 5
# for stress controlled defined CSR > 0
CSR = 0.1
# for strain controlled defined Epslim > 0
Epslim = 0.
# HSS parameters
t.model['stress-1-effective']= -sig3/k0
t.model['stress-2-effective']= -sig3
t.model['stress-3-effective']= -sig3

t.CyDSS(
    sig3, 
    drainage,
    Ncycles=Ncycles, 
    Epslim=Epslim, 
    CSR=CSR,
    figname = r'.\output\CyDSS_res.png', 
    filename= r'.\output\CyDSS_res.csv')
''' 