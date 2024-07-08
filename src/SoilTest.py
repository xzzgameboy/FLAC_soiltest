'''
Copyright SRK Consulting (Canada) Inc. 2022 - All Rights Reserved
Licensed under CC BY-SA https://creativecommons.org/licenses/by-sa/4.0/
Author: info@srk.com
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import csv

import itasca as it

class SoilTest:
        
    def __init__(self):
        
        self.model = None
        self.thread = 1
        
        # fluid properties
        self.kw = 2.e6
        self.u0 = 0.
        
        # numerical parameters
        self.strain_step = 1.e-4
        self.sub_steps = 50
        self.strain_rate = self.strain_step/self.sub_steps
        self.tol = 1.e-4
        
        # control points
        self.z = None
        self.gp = None
        
        # tracking variables
        self.epsxx = None
        self.epsyy = None
        self.epszz = None
        self.epsxy = None
        self.epsxz = None
        self.epsyz = None
        
        self.sigxx = None
        self.sigyy = None
        self.sigzz = None
        self.sigxy = None
        self.sigxz = None
        self.sigyz = None

        self.pwater = None

        self.step = None
        
    def set_up (self):
    
        it.command('python-reset-state false')
        
        new_template = ('''
            project new
            model deterministic on
            fish automatic-create off
            model large-strain off
            program threads {thread}
            ''')
        new_command = new_template.format(thread = self.thread)
            
        it.command(new_command)
        
        it.command('''
            model config fluid
            model fluid active off
            ''')
            
        it.command('zone create brick size 1 1 1')
        
        # control points
        center = (0.5, 0.5, 0.5)
        corner = (0., 0., 0.)

        self.z = it.zone.near(center)
        self.gp = it.gridpoint.near(corner)
        
    def stress_ini (self, sig3, k0=1.):
    
        sigxx = sig3
        sigyy = sig3
        sigzz = sig3/k0
        
        # initial stress state
        ini_template = '''
        zone initialize stress xx -{sigxx} yy -{sigyy} zz -{sigzz}
        '''
        ini_assign = ini_template.format(sigxx=sigxx, sigyy=sigyy, sigzz=sigzz)
        it.command(ini_assign)
    
    
    def set_drainage (self, drainage):
    
        if drainage == 'undrained':
        
            kw = self.kw
            
        elif drainage == 'drained':  
        
            kw = 0.
            
        fluid_template = '''
            zone fluid cmodel assign {model}
            zone gridpoint initialize fluid-modulus {kw}
            zone gridpoint initialize pore-pressure {u0}
            '''
        fluid_asign = fluid_template.format(model='isotropic', 
                kw=kw, u0=self.u0)
        it.command(fluid_asign)
 
           
    def assign_model (self):
    
        # constitutive model
        for key in self.model:
            if key == 'model':
                it.command('zone cmodel assign {}'.format(self.model[key]))
            else:
                it.command('zone property {} {}'.format(key, self.model[key]))
                
    def hist_ini (self):
        
        stress = self.z.stress()    
        self.sigxx = np.array([stress.xx()])
        self.sigyy = np.array([stress.yy()])
        self.sigzz = np.array([stress.zz()])
        self.sigxy = np.array([stress.xy()])
        self.sigxz = np.array([stress.xz()])
        self.sigyz = np.array([stress.yz()])
        
        strain = self.z.strain()
        self.epsxx = np.array([strain.xx()])
        self.epsyy = np.array([strain.yy()])
        self.epszz = np.array([strain.zz()])
        self.epsxy = np.array([strain.xy()])
        self.epsxz = np.array([strain.xz()])
        self.epsyz = np.array([strain.yz()])
        
        self.step = np.array([it.cycle()])
        
        self.pwater = np.array([self.gp.pp()])    
        
    def hist_update (self):
    
        stress = self.z.stress()      
        
        self.sigxx = np.append(self.sigxx, stress.xx())
        self.sigyy = np.append(self.sigyy, stress.yy())
        self.sigzz = np.append(self.sigzz, stress.zz())
        self.sigxy = np.append(self.sigxy, stress.xy())
        self.sigxz = np.append(self.sigxz, stress.xz())
        self.sigyz = np.append(self.sigyz, stress.yz())
        
        strain = self.z.strain()
        self.epsxx = np.append(self.epsxx, strain.xx())
        self.epsyy = np.append(self.epsyy, strain.yy())
        self.epszz = np.append(self.epszz, strain.zz())
        self.epsxy = np.append(self.epsxy, strain.xy())
        self.epsxz = np.append(self.epsxz, strain.xz())
        self.epsyz = np.append(self.epsyz, strain.yz())
        
        self.step = np.append(self.step, it.cycle())
        
        self.pwater = np.append(self.pwater, self.gp.pp())
        
    def calc_TXX (self, epsa_max=5.):
                      
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0 = stress.xx(), stress.yy()
        strain = self.z.strain()
        epsa0 = strain.zz()
        
        # boundary conditions             
        bc_template =''' 
        zone gridpoint fix velocity-z
        zone face apply stress-xx {sigxx0} range union position-x 0 position-x 1
        zone face apply stress-yy {sigyy0} range union position-y 0 position-y 1
        '''
        bc_assign = bc_template.format(sigxx0=sigxx0, sigyy0=sigyy0)
        it.command(bc_assign)
        
        # run model
        run_template = '''
        zone face apply velocity-z -{strain_rate} range position-z 1.0
        model step {step}
        zone face apply velocity-z 0.0 range position-z 1.0
        model solve ratio {tol}
        '''
        run_assign = run_template.format(
            step=self.sub_steps, 
            strain_rate=self.strain_rate, 
            tol=self.tol)

        eps_a = epsa0
        while eps_a <= epsa_max:
            it.command(run_assign)
            self.hist_update()
            strain = self.z.strain()
            eps_a = -strain.zz()*100.
        
    def calc_TXXcv (self, epsa_max=5.):
    
        self.set_drainage ('drained')
                      
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        p0 = (sigxx0 + sigyy0 + sigzz0)/3.
        u0 = self.gp.pp()
        q0 = sigzz0 - sigxx0
        
        strain = self.z.strain()
        epsa0 = strain.zz()
        
        # boundary conditions             
        bc_template =''' 
        zone gridpoint fix velocity-z
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-x
        '''
        bc_assign = bc_template.format(sigxx0=sigxx0, sigyy0=sigyy0)
        it.command(bc_assign)
        
        run_template = '''
        zone face apply velocity-z -{strain_rate} range position-z 1.0
        zone face apply velocity-y {half_sr} range position-y 1.0
        zone face apply velocity-x {half_sr} range position-x 1.0
        model step {step}
        zone face apply velocity-z 0.0 range position-z 1.0
        zone face apply velocity-y 0.0 range position-y 1.0
        zone face apply velocity-x 0.0 range position-x 1.0
        model solve ratio {tol}
        '''
        run_assign = run_template.format(
            step = self.sub_steps, 
            strain_rate = self.strain_rate, 
            half_sr = self.strain_rate/2.,
            tol = self.tol)
        
        eps_a = epsa0
        while eps_a <= epsa_max:
            it.command(run_assign)
            self.hist_update()
            # calculate pore pressure and total stresses
            stress = self.z.stress()
            q = stress.zz() - stress.xx()
            p = (stress.xx() + stress.yy() + stress.zz())/3.
            delta_q = q - q0
            p_tot = p0 + delta_q/3.
            delta_u = -(p_tot - p)
            u = u0 + delta_u
            self.pwater[-1] = u
            self.sigxx[-1] += -u
            self.sigyy[-1] += -u
            self.sigzz[-1] += -u
            # update axial strain
            strain = self.z.strain()
            eps_a = -strain.zz()*100.
            
            
    def TXX (self, sig3, k0, drainage, epsa_max,
        savefig=True, savefile=True,
        figname = r'.\output\TXX_res.png', filename=r'.\output\TXX_res.csv'):
        
        # calc triaxial
        self.set_up()
        self.assign_model()
        self.stress_ini(sig3=sig3, k0=k0)
        self.hist_ini()
        if drainage == 'constant volume':
            self.calc_TXXcv(epsa_max=epsa_max)
        else:
            self.set_drainage(drainage)
            self.calc_TXX(epsa_max=epsa_max)
        
        # plot results
        if savefig:
            self.plot_ax(zero_origin = True, figname=figname)
        
        # save results
        if savefile:
            self.save_res (filename=filename)
        
            
    def calc_OED (self, load):
        
        # boundary conditions             
        it.command(''' 
        zone face apply velocity-y 0. range union position-y 0 position-y 1
        zone face apply velocity-x 0. range union position-x 0 position-x 1
        zone face apply velocity-z 0. range position-z 0
        ''')
     
        # run model
        run_template = '''
        zone face apply stress-zz -{load} range position-z 1.0
        model solve ratio {tol}
        '''
        
        for l in load:
            run_assign = run_template.format(
                load=l, 
                tol=self.tol)
            it.command(run_assign)
            self.hist_update()
            
    def OED (self, load, sig0=0.,
        savefig=True, savefile=True,
        figname = r'.\output\OED_res.png', filename=r'.\output\OED_res.csv'):
    
        self.set_up()
        self.assign_model()
        self.stress_ini(sig3=sig0)
        self.hist_ini()
        self.set_drainage('drained')
        self.calc_OED(load)
        
        # plot results
        plt.plot(-self.sigzz, self.epszz*100,  marker='o')
        plt.xscale('log')
        plt.grid(visible=True, which='both')
        plt.xlabel('Vertical effective stress (kPa)')
        plt.ylabel('Volumetric strain (%)')
        plt.tight_layout()

        if savefig:
            plt.savefig(figname)
        
        # save results
        if savefile:
            self.save_res (filename=filename)


    def calc_CyTXX (self, Ncycles = 5, Epslim = 0.0, CSR = 0.0):
        
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        p_ini = np.absolute( (sigxx0 + sigyy0 + sigzz0)/3 )
        strain = self.z.strain()
        epsa0 = strain.zz()  
        
        #stress limit and strain limit
        if ((Epslim > 0) and (CSR==0)):
        
            mode = 'StrainControl'
            var_limit = Epslim
            
        elif ((Epslim == 0) and (CSR >0)):
        
            mode = 'StressControl'
            var_limit = CSR * p_ini * 2.
            
        else:
        
            mode = 'Error'
        
        # boundary conditions
        bc_template =''' 
        zone gridpoint fix velocity-z range position-z -0.2 0.2
        zone face apply stress-xx {sigxx0} range union position-x 0 position-x 1
        zone face apply stress-yy {sigyy0} range union position-y 0 position-y 1
        zone face apply stress-zz {sigzz0} range position-z 0.8 1.2'''
        bc_assign = bc_template.format(sigxx0=sigxx0, sigyy0=sigyy0, 
            sigzz0=sigzz0)
        it.command(bc_assign)

        # run model
        run_template = '''
        zone face apply velocity-z {strain_rate} range position-z 1.0
        model step {step}
        zone face apply velocity-z 0.0 range position-z 1.0
        model solve ratio {tol}
        '''
        
        theta = 0.
        dtheta= theta
        Ttheta= 0.
        SignLoad = -1.0
        eps_a = epsa0       

        if mode != 'Error':
        
            while Ttheta < 2 * np.pi * Ncycles:

                run_assign = run_template.format(
                step = self.sub_steps, 
                strain_rate = self.strain_rate * SignLoad, 
                tol = self.tol)

                it.command(run_assign)
                
                ## reading and saving stress and stress values
                self.hist_update()    
                            
                if mode == 'StrainControl':
                
                    strain = self.z.strain()
                    epszz = strain.zz()       
                    var = -epszz*100.

                elif mode == 'StressControl':
                
                    stress = self.z.stress()
                    sigxx, sigzz = stress.xx(), stress.zz()
                    var = sigzz - sigxx
                   
                dtheta = theta
                
                if np.absolute(var) > var_limit:

                    SignLoad = SignLoad * -1.
                    theta = np.arcsin(var/np.absolute(var))

                else:
                
                    theta = np.arcsin(var/np.absolute(var_limit))
                 
                Ttheta = Ttheta + np.absolute(theta - dtheta)
                
        else:
        
            print('Error, review inputs')
            
    def calc_CyTXXcv(self, Ncycles = 5, Epslim = 0.0, CSR = 0.0):
        
        self.set_drainage ('drained')
                      
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        p0 = (sigxx0 + sigyy0 + sigzz0)/3.
        u0 = self.gp.pp()
        q0 = sigzz0 - sigxx0
        
        strain = self.z.strain()
        epsa0 = strain.zz()
        
        #stress limit and strain limit
        if ((Epslim > 0) and (CSR==0)):
        
            mode = 'StrainControl'
            var_limit = Epslim
            
        elif ((Epslim == 0) and (CSR >0)):
        
            mode = 'StressControl'
            var_limit = CSR * p0 * 2.
            
        else:
        
            mode = 'Error'
        
        # boundary conditions
        bc_template =''' 
        zone gridpoint fix velocity-z
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-x
        '''
        bc_assign = bc_template.format(sigxx0=sigxx0, sigyy0=sigyy0, 
            sigzz0=sigzz0)
        it.command(bc_assign)

        # run model
        run_template = '''
        zone face apply velocity-z {strain_rate} range position-z 1.0
        zone face apply velocity-y {half_sr} range position-y 1.0
        zone face apply velocity-x {half_sr} range position-x 1.0
        model step {step}
        zone face apply velocity-z 0.0 range position-z 1.0
        zone face apply velocity-y 0.0 range position-y 1.0
        zone face apply velocity-x 0.0 range position-x 1.0
        model solve ratio {tol}
        '''
        
        theta = 0.
        dtheta= theta
        Ttheta= 0.
        SignLoad = -1.0
        eps_a = epsa0       

        if mode != 'Error':
        
            while Ttheta < 2 * np.pi * Ncycles:

                run_assign = run_template.format(
                step = self.sub_steps, 
                strain_rate = self.strain_rate * SignLoad,
                half_sr = -1*SignLoad*self.strain_rate/2.,
                tol = self.tol)

                it.command(run_assign)
                
                ## reading and saving stress and stress values
                self.hist_update()    
                
                # calculate pore pressure and total stresses
                stress = self.z.stress()
                q = stress.zz() - stress.xx()
                p = (stress.xx() + stress.yy() + stress.zz())/3.
                delta_q = q - q0
                p_tot = p0 + delta_q/3.
                delta_u = -(p_tot - p)
                u = u0 + delta_u
                self.pwater[-1] = u
                self.sigxx[-1] += -u
                self.sigyy[-1] += -u
                self.sigzz[-1] += -u
                
                if mode == 'StrainControl':
                
                    strain = self.z.strain()
                    epszz = strain.zz()       
                    var = -epszz*100.

                elif mode == 'StressControl':
                
                    var = q
                   
                dtheta = theta
                
                if np.absolute(var) > np.absolute(var_limit):
                    SignLoad = SignLoad * -1.
                    theta = np.arcsin(var/np.absolute(var))

                else:
                    theta = np.arcsin(var/np.absolute(var_limit))
                 
                Ttheta = Ttheta + np.absolute(theta - dtheta)
                
        else:
        
            print('Error, review inputs')
            
    def CyTXX (self, sig3, drainage, 
        k0=1.0, Ncycles = 5, Epslim = 0.0, CSR = 0.0,
        savefig=True, savefile=True,
        figname = r'.\output\CyTXX_res.png', 
        filename=r'.\output\CyTXX_res.csv'):
    
        self.set_up()
        self.assign_model()
        self.stress_ini(sig3=sig3, k0=k0)
        self.hist_ini()
        
        
        if drainage == 'constant volume':
            self.calc_CyTXXcv(Ncycles=Ncycles, Epslim=Epslim, CSR=CSR)
        else:
            self.set_drainage(drainage)
            self.calc_CyTXX(Ncycles=Ncycles, Epslim=Epslim, CSR=CSR)
        
        # plot results
        if savefig:
            self.plot_ax(zero_origin = False, figname=figname)
        
        # save results
        if savefile:
            self.save_res (filename=filename)
            
    def calc_DSS (self, epsyz_max=5.):
                      
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        strain = self.z.strain()
        epsyz0 = strain.yz()
        
        # boundary conditions               
        bc_template =''' 
        zone gridpoint fix velocity-x
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-z range position-z 0
        zone attach gridpointid 4 to-gridpointid 8 snap off
        zone attach gridpointid 6 to-gridpointid 8 snap off
        zone attach gridpointid 7 to-gridpointid 8 snap off
        zone face apply stress-zz {sigzz0} range position-z 1
        '''
        bc_assign = bc_template.format(sigzz0=sigzz0)
        it.command(bc_assign)
        
        # run model
        run_template = '''
        zone face apply velocity-y -{strain_rate} range position-z 1
        model step {step}
        zone face apply velocity-y 0.0 range position-z 1
        model solve ratio {tol}
        '''
        run_assign = run_template.format(
            step=self.sub_steps, 
            strain_rate=self.strain_rate, 
            tol=self.tol)

        epsyz = epsyz0*100.
        while abs(epsyz) <= epsyz_max:
            it.command(run_assign)
            self.hist_update()
            strain = self.z.strain()
            epsyz = strain.yz()*100.
    
    def calc_DSScv (self, epsyz_max=5.):
        
        self.set_drainage ('drained')
        
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        strain = self.z.strain()
        epsyz0 = strain.yz()
        u0 = self.gp.pp()
        
        # boundary conditions         
        bc_assign =''' 
        zone gridpoint fix velocity-z
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-x
        '''
        it.command(bc_assign)
        
        # run model
        run_template = '''
        zone face apply velocity-y -{strain_rate} range position-z 1
        model step {step}
        zone face apply velocity-y 0.0 range position-z 1
        model solve ratio {tol}
        '''
        run_assign = run_template.format(
            step=self.sub_steps, 
            strain_rate=self.strain_rate, 
            tol=self.tol)

        epsyz = epsyz0*100.
        while abs(epsyz) <= epsyz_max:
            it.command(run_assign)
            self.hist_update()
            # calculate pore pressure and total stresses
            stress = self.z.stress()
            delta_u = stress.zz() - sigzz0
            u = u0 + delta_u
            self.pwater[-1] = u
            self.sigxx[-1] += -u
            self.sigyy[-1] += -u
            self.sigzz[-1] += -u
            # upsate shear strain
            strain = self.z.strain()
            epsyz = strain.yz()*100.
            
    def DSS (self, sig3, k0, drainage, epsyz_max,
        savefig=True, savefile=True,
        figname = r'.\output\DSS_res.png', filename=r'.\output\DSS_res.csv'):
        
        self.set_up()
        self.assign_model()
        self.stress_ini(sig3=sig3, k0=k0)
        self.hist_ini()
        if drainage == 'constant volume':
            self.calc_DSScv(epsyz_max)
        else:
            self.set_drainage(drainage)
            self.calc_DSS(epsyz_max)
        
         # plot results
        if savefig:
            self.plot_ps(figname=figname)
        
        # save results
        if savefile:
            self.save_res (filename=filename)
        
    def calc_CyDSS(self, NCycles = 5, Epslim = 0.0, CSR = 0.0):
                      
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        strain = self.z.strain()
        epsyz0 = strain.yz()
 
        #stress limit and strain limit
        if ((Epslim > 0) and (CSR==0)):
            mode = 'StrainControl'
            var_limit = Epslim
        elif ((Epslim == 0) and (CSR >0)):
            mode = 'StressControl'
            var_limit = CSR * sigzz0 
        else:
            mode = 'Error'
 
        # boundary conditions         
        bc_template = ''' 
        zone gridpoint fix velocity-x
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-z range position-z 0
        zone attach gridpointid 4 to-gridpointid 8 snap off
        zone attach gridpointid 6 to-gridpointid 8 snap off
        zone attach gridpointid 7 to-gridpointid 8 snap off
        zone face apply stress-zz {sigzz0} range position-z 1
        '''
        bc_assign = bc_template.format(sigzz0=sigzz0)
        it.command(bc_assign)
        
        # run model
        run_template = '''
        zone face apply velocity-y {strain_rate} range position-z 1
        model step {step}
        zone face apply velocity-y 0.0 range position-z 1
        model solve ratio {tol}
        '''
        
        theta = 0.
        dtheta= theta
        Ttheta= 0.
        SignLoad = -1.0        
        epsyz = epsyz0*100.
        
        if mode != 'Error':
        
            while Ttheta < 2 * np.pi * NCycles:
                
                run_assign = run_template.format(
                step=self.sub_steps, 
                strain_rate=self.strain_rate * SignLoad, 
                tol=self.tol)                
                it.command(run_assign)
                
                ## reading and saving stress and stress values
                self.hist_update()    
                
                if mode == 'StrainControl':
                    strain = self.z.strain()
                    epsyz = strain.yz()       
                    var = -epsyz*100.
                elif mode == 'StressControl':
                    stress = self.z.stress()
                    sigyz= stress.yz()
                    var = sigyz
                
                dtheta = theta
                
                if np.absolute(var) > np.absolute(var_limit):
                    SignLoad = SignLoad * -1.
                    theta = np.arcsin(var/np.absolute(var))
                else:
                    theta = np.arcsin(var/np.absolute(var_limit))
                 
                Ttheta = Ttheta + np.absolute(theta - dtheta)
        else:
            print('Error, review inputs')
            
    def calc_CyDSScv(self, NCycles = 5, Epslim = 0.0, CSR = 0.0):
        
        self.set_drainage ('drained')
        
        # read initial stress state
        stress = self.z.stress()
        sigxx0, sigyy0, sigzz0 = stress.xx(), stress.yy(), stress.zz()
        strain = self.z.strain()
        epsyz0 = strain.yz()
        p0 = (sigxx0 + sigyy0 + sigzz0)/3.
        u0 = self.gp.pp()
        q0 = sigzz0 - sigxx0
        
        #stress limit and strain limit
        if ((Epslim > 0) and (CSR==0)):    
            mode = 'StrainControl'
            var_limit = Epslim
        elif ((Epslim == 0) and (CSR >0)):       
            mode = 'StressControl'
            var_limit = CSR * sigzz0           
        else:      
            mode = 'Error'

        # boundary conditions         
        bc_assign =''' 
        zone gridpoint fix velocity-z
        zone gridpoint fix velocity-y
        zone gridpoint fix velocity-x
        '''
        it.command(bc_assign)
        
        # run model
        run_template = '''
        zone face apply velocity-y {strain_rate} range position-z 1
        model step {step}
        ;zone face apply velocity-y 0.0 range position-z 1
        model solve ratio {tol}
        '''       
        
        theta = 0.
        dtheta= theta
        Ttheta= 0.
        SignLoad = -1.0        
        epsyz = epsyz0*100.
        
        if mode != 'Error':
        
            while Ttheta < 2 * np.pi * NCycles:
                
                run_assign = run_template.format(
                step= self.sub_steps, 
                strain_rate= self.strain_rate * SignLoad,
                tol= self.tol)                
                it.command(run_assign)
                
                ## reading and saving stress and stress values
                self.hist_update()    
                
                # calculate pore pressure and total stresses
                stress = self.z.stress()
                delta_u = stress.zz() - sigzz0
                u = u0 + delta_u
                self.pwater[-1] = u
                self.sigxx[-1] += -u
                self.sigyy[-1] += -u
                self.sigzz[-1] += -u
                
                if mode == 'StrainControl':               
                    strain = self.z.strain()
                    epsyz = strain.yz()       
                    var = epsyz*100.
                elif mode == 'StressControl':               
                    stress = self.z.stress()
                    sigyz= stress.yz()
                    var = sigyz                    
                
                dtheta = theta

                if np.absolute(var) > np.absolute(var_limit):
                    SignLoad = SignLoad * -1.
                    theta = np.arcsin(var/np.absolute(var))                   
                else:                   
                    theta = np.arcsin(var/np.absolute(var_limit))
                 
                Ttheta = Ttheta + np.absolute(theta - dtheta)                
        else:       
            print('Error, review inputs')
            
    def CyDSS(self, sig3, drainage, 
        k0=1.0, Ncycles = 5, Epslim = 0.0, CSR = 0.0,
        savefig=True, savefile=True,
        figname = r'.\output\CyDSS_res.png', 
        filename=r'.\output\CyDSS_res.csv'):
    
        self.set_up()
        self.assign_model()
        self.stress_ini(sig3=sig3, k0=k0)
        self.hist_ini()
        
        if drainage == 'constant volume':
            self.calc_CyDSScv(NCycles=Ncycles, Epslim=Epslim, CSR=CSR)
        else:
            self.set_drainage(drainage)
            self.calc_CyDSS(NCycles=Ncycles, Epslim=Epslim, CSR=CSR)
        
        # plot results
        if savefig:
            self.plot_ps(zero_origin = False, figname=figname)
        
        # save results
        if savefile:
            self.save_res (filename=filename)
            
    def plot_ax (self, zero_origin = True, figname = r'.\output\res.png'):
    
        # calculate key strain and stress variables
        u = self.pwater
        q = -(self.sigzz - self.sigxx)
        p = -(self.sigzz + self.sigxx + self.sigyy)/3. - u
        epsv = -(self.epszz + self.epsxx + self.epsyy)*100.
        epsa = -self.epszz*100.

        fig_size = (10, 10)
        fig = plt.figure(figsize = fig_size)
        gs = GridSpec(2, 2)
        
        # plot 01 - q vs epsa 
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(epsa, q)
        if zero_origin: ax1.set_ylim([0, None])
        ax1.set_xlabel('Axial Strain, $\epsilon_a$ (%)')
        ax1.set_ylabel('Deviatoric Stress, q (kPa)')
        ax1.grid(False)
        
        # plot 02 - q vs p
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(p, q)
        if zero_origin: 
            ax2.set_ylim([0, None])
            ax2.set_xlim([0, None])
        ax2.set_xlabel('Mean effective stress,  p (kPa)')
        ax2.set_ylabel('Deviatoric Stress, q (kPa)')
        ax2.grid(False)

        # plot - epsv vs epsa 
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.plot(epsa, epsv)
        ax3.set_xlabel('Axial strain, $\epsilon_a$ (%)')
        ax3.set_ylabel('Volumetric strain, $\epsilon_v$ (%)')
        
        # plot - u vs epsa
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.plot(epsa, u)
        ax4.set_xlabel('Axial strain, $\epsilon_a$ (%)')
        ax4.set_ylabel('Pore pressure, $p_w$ (kPa)')
        
        # save plot
        fig.savefig(figname, bbox_inches='tight')
        
    def plot_ps (self, zero_origin = True, figname = r'.\output\res.png'):
    
        # calculate key strain and stress variables
        u = self.pwater
        tau = -self.sigyz
        sigv = -self.sigzz - u
        epsyz = -self.epsyz*100.
        epsv = -(self.epszz + self.epsxx + self.epsyy)*100.

        fig_size = (10, 10)
        fig = plt.figure(figsize = fig_size)
        gs = GridSpec(2, 2)
        
        # plot 01 - q vs epsa 
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(epsyz, tau)
        if zero_origin: ax1.set_ylim([0, None])
        ax1.set_xlabel('Shear strain, $\epsilon_{yz}$ (%)')
        ax1.set_ylabel(r'Shear Stress, $\tau$ (kPa)')
        ax1.grid(False)
        
        # plot 02 - q vs p
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(sigv, tau)
        if zero_origin: 
            ax2.set_ylim([0, None])
            ax2.set_xlim([0, None])
        ax2.set_xlabel('Vertical effective stress,  $\sigma_v$ (kPa)')
        ax2.set_ylabel(r'Shear Stress, $\tau$ (kPa)')
        ax2.grid(False)

        # plot - epsv vs epsa 
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.plot(epsyz, epsv)
        ax3.set_xlabel('Shear strain, $\epsilon_{yz}$ (%)')
        ax3.set_ylabel('Volumetric strain, $\epsilon_v$ (%)')
        
        # plot - u vs epsa
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.plot(epsyz, u)
        ax4.set_xlabel('Shear strain, $\epsilon_{yz}$ (%)')
        ax4.set_ylabel('Pore pressure, $p_w$ (kPa)')
        
        # save plot
        fig.savefig(figname, bbox_inches='tight')
    
    def save_res (self, filename = r'.\output\res.csv'):
   
        with open(filename,"w") as file:
    
            #write first line (header) in the file
            file.writelines((
                '{}, ' +
                '{}, {}, {}, {}, {}, {}, '+
                '{}, {}, {}, {}, {}, {}, '+
                '{}\n').
                format(
                'Step', 
                'sigx', 'sigy', 'sigz', 'sigxy', 'sigxz', 'sigyz',
                'epsx', 'epsy', 'epsz', 'epsxy', 'epsxz', 'epsyz',
                'pore pressure'))
            
            #write results
            file.writelines([(
                '{}, ' +
                '{}, {}, {}, {}, {}, {}, '+
                '{}, {}, {}, {}, {}, {}, '+
                '{}\n').
                format(
                st, 
                sx, sy, sz, sxy, sxz, syz,
                stx, sty, stz, stxy, stxz, styz,
                pwater)
                #loop to write results
                for (
                st, 
                sx, sy, sz, sxy, sxz, syz,
                stx, sty, stz, stxy, stxz, styz,
                pwater) in
                zip(
                self.step,
                self.sigxx, self.sigyy, self.sigzz, self.sigxy, self.sigxz, 
                self.sigyz,
                self.epsxx, self.epsyy, self.epszz, self.epsxy, self.epsxz, 
                self.epsyz,
                self.pwater)
            ])
            
    def read_labdata (self, path):
        
        #reading .csv file
        with open(path, "r") as f:
            csv_reader = csv.DictReader(f)
            data = list(csv_reader)

        #collecting all data in one variable
        test_ID_all = [data[i]['ï»¿test_ID'] for i in range(len(data))]
        eps_a_all   = [float(data[i]['strain_axial']) for i in range(len(data))]
        q_all       = [float(data[i]['q_kPa']) for i in range(len(data))]
        p_all       = [float(data[i]['p_kPa']) for i in range(len(data))]
        eps_v_all   = [float(data[i]['strain_vol']) for i in range(len(data))]
        u_all       = [float(data[i]['u_kPa']) for i in range(len(data))]

        #obtain size within each variable
        n = 0
        j = 0
        ini = test_ID_all[0]
        id = [ini]
        index = []
        for ID in test_ID_all :
            j = j+1
            if ID == ini :
                n = n+1
            else:
                id.append(ID)
                index.append(n)
                ini = ID
                n= 1
            if j == len(test_ID_all) :
                index.append(n)
                
        #index position
        pos = np.cumsum(index)

        #obtain variables for each test
        eps_a = []
        for i in range(len(pos)) :
            if i==0 :
                eps_a.append(eps_a_all[0:pos[i]])
            else:
                eps_a.append(eps_a_all[pos[i-1]:pos[i]])

        q = []
        for i in range(len(pos)) :
            if i==0 :
                q.append(q_all[0:pos[i]])
            else:
                q.append(q_all[pos[i-1]:pos[i]])

        p = []
        for i in range(len(pos)) :
            if i==0 :
                p.append(p_all[0:pos[i]])
            else:
                p.append(p_all[pos[i-1]:pos[i]])
                
        eps_v = []
        for i in range(len(pos)) :
            if i==0 :
                eps_v.append(eps_v_all[0:pos[i]])
            else:
                eps_v.append(eps_v_all[pos[i-1]:pos[i]])
    
        u = []
        for i in range(len(pos)) :
            if i==0 :
                u.append(u_all[0:pos[i]])
            else:
                u.append(u_all[pos[i-1]:pos[i]])
        
        eps_a = [np.array(x)*100 for x in eps_a]
        eps_v = [np.array(x)*100 for x in eps_v]
        
        return eps_a, q, p, eps_v, u
