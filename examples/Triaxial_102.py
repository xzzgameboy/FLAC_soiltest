'''
Copyright SRK Consulting (Canada) Inc. 2022 - All Rights Reserved
Licensed under CC BY-SA https://creativecommons.org/licenses/by-sa/4.0/
Author: info@srk.com

================================================================================
Fundao Triaxial Test

    - Runs isotropic triaxial tests from Fundao's failure
    - Uses Norsand
    
================================================================================
'''

from src import SoilTest as st
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

t = st.SoilTest()

# ======================================================================== INPUT

tests = [
    ['TX-01', 'drained', 0.808, 300.],
    ['TX-02', 'undrained', 0.838, 200.],
    ['TX-04', 'undrained', 0.796, 600.],
    ['TX-06', 'undrained', 0.804, 300.],
    ['TX-08', 'drained', 0.728, 300.],
    ['TX-10', 'drained', 0.595, 600.],
    ['TX-11', 'drained', 0.611, 200.],
    ['TX-12', 'drained', 0.69, 400.],
    ['TX-13', 'undrained', 0.689, 600.],
    ['TX-18', 'drained', 0.809, 100.],
    ['TX-19', 'drained', 0.882, 50.],
    ['TX-29', 'drained', 0.741, 400.],
    ['TX-30', 'drained', 0.784, 200.]
    ]
    
    
# maximum axial strain in percentage
max_strain = 30.

# material parameters
t.model = {
    'model'             : 'norsand',
    # basic properties 
    'critical-state-1'  : 0.80,# Gamma when C3 = 0
    'critical-state-2'  : 0.03,# lambda_e when C3 = 0
    'critical-state-3'  : 0.65,# 0 by default
    'factor-coupling'   : 0.38, # N, volumetric coupling factor
    'factor-dilatancy'  : 7.3, # \chi_tc
    'hardening-0'       : 156., # H0, plastic hardening modulus when psi=0
    'shear-reference'   : 20.e3,# Gref, reference shear modulus
    'ratio-critical'    : 1.33,# Mtc, stress ratio at the critical state line 
    # Advanced properties
    'exponent'          : 0.4, # m, pressure exponent for elasticity
    'over-consolidation-ratio': 1.0, # OCR
    'hardening-y'       : 756., # H_psi, hardening law dependance on psi
    'pressure-reference': 100., # pref
    'poisson'           : 0.3, # nu
    # Initial properties
    'stress-xy-initial' : 0., 
    'stress-yz-initial' : 0.,
    'stress-xz-initial' : 0.
    } 
 
# reading lab data
path = 'examples/data/test_results.csv'

lab_epsa, lab_q, lab_p, lab_epsv, lab_u = t.read_labdata(path)

# ===================================================================== run test

for i in range(len(tests)):

    start = time.time()
    
    # define testID, drainage, void ratio and sig_mean for each test
    testID, drainage, void, sig_mean = tests[i]
    
    # assign initial conditions
    t.model['void-initial'] = void
    t.model['stress-xx-initial'] = -sig_mean
    t.model['stress-yy-initial'] = -sig_mean
    t.model['stress-zz-initial'] = -sig_mean
    
    # activates softening flag if undrained
    if drainage == 'drained':
    
        t.model['index-softening'] = 0
        
    elif drainage == 'undrained':
    
        t.model['index-softening'] = 1
    
    t.TXX(
        sig3 = sig_mean, 
        k0 = 1.0, 
        drainage = drainage, 
        epsa_max = max_strain,
        savefig = False, 
        savefile = False,)

    # ============================================================= plot results
    fig_size = (15, 5)
    fig = plt.figure(figsize = fig_size)
    gs = GridSpec(1, 3)

    u = t.pwater
    q = -(t.sigzz - t.sigxx)
    p = -(t.sigzz + t.sigxx + t.sigyy)/3. - u
    epsv = -(t.epszz + t.epsxx + t.epsyy)*100.
    epsa = -t.epszz*100.

    # plot 01 - q vs epsa 
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(epsa, q, label='FLAC3D')
    ax1.plot(lab_epsa[i], lab_q[i], label='Test Results')
    ax1.set_xlabel('Axial Strain, $\epsilon_a$ (%)')
    ax1.set_ylabel('Deviatoric Stress, q (kPa)')
    ax1.set_xlim([0, max_strain])
    ax1.set_ylim([0, None])
    ax1.grid(False)
    ax1.legend()

    # plot 02 - q vs p
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(p, q, label='FLAC3D')
    ax2.plot(lab_p[i], lab_q[i], label='Test Results')
    ax2.set_xlabel('Mean effective stress,  p (kPa)')
    ax2.set_ylabel('Deviatoric Stress, q (kPa)')
    ax2.set_xlim([0, None])
    ax2.set_ylim([0, None])
    ax2.grid(False)
    ax2.legend()

    # plot 03 - undrained
    ax3 = fig.add_subplot(gs[0, 2])
    if drainage=='drained':
        # plot - epsv vs epsa 
        ax3.plot(epsa, epsv*-1, label='FLAC3D')
        ax3.plot(lab_epsa[i], lab_epsv[i], label='Test Results')
        ax3.set_xlabel('Axial strain, $\epsilon_a$ (%)')
        ax3.set_ylabel('Volumetric strain, $\epsilon_v$ (%)')
        ax3.set_xlim([0, max_strain])
    # plot 02 - drained    
    elif drainage== 'undrained':
        # plot - u vs epsa
        ax3.plot(epsa, u, label='FLAC3D')
        ax3.plot(lab_epsa[i], lab_u[i], label='Test Results')
        ax3.set_xlabel('Axial strain, $\epsilon_a$ (%)')
        ax3.set_ylabel('Pore pressure, $p_w$ (kPa)')
        ax3.set_xlim(0, None)
    
    ax3.grid(False)
    ax3.legend()
    
    # title
    title = testID +'; '+drainage+'; void ratio: '+str(void)+'; confining stress: '+str(sig_mean)+' kPa'
    plt.suptitle(title, fontsize=14, y=1.05)
    fig.tight_layout()
    
    # fig.title(testID)
    figname = r'.\output\Triaxial_102_'+testID+'.pdf'
    fig.savefig(figname, bbox_inches='tight')

    # ============================================================= save results
    filename = r'.\output\Triaxial_102_'+testID+'.csv'
    t.save_res(filename)