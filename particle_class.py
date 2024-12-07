import numpy as np
import math




e = 1.6021766e-19
eta = 1e-3 #Pa
epsilon0 = 8.85541878e-12
epsilon_r = 80.10
epsilon = epsilon_r*epsilon0
Na = 6.022e23
kb = 1.3806488e-23

class _particle_():
    def __init__(self, T, m,c, index):
        I = c*Na
        kappa = 2.2912074e-3*np.sqrt(I/T)
        phi0 = 0.1
        self.R = 1e-5 #m
        self.kappa = kappa
#bulk DLVO

#modified DLVO
        z = np.tanh(e*phi0/(4*kb*T))
        self.DLVO = 64*I*kb*T*(z**2)/(self.R)
#Van der Waals
        self.VanderWaals = 0.95e-20*self.kappa**3/(12*np.pi)
        self.T = T
#assumin partícules esfèriques
        self.D = kb*self.T/(6*np.pi*eta*self.R)
        density = 2.12e3
        self.m = 4/3*np.pi*(self.R)**3*density 
        self.index = index
        self.lengh = 2*np.pi/(2.64)*5e-3*kappa #wavelenght 
        sigma = np.sqrt(kb*T/m)
        print('sigma', sigma)
        self.rdot = (np.array([np.random.normal(loc=0, scale = sigma), np.random.normal(loc=0, scale = sigma)]))/(2*self.D*self.kappa)
        print(self.rdot)
#        self.pos = np.array([np.random.randint(low = 0, high = 2)*self.lengh+np.random.normal(loc=0, scale = self.lengh/20), np.random.normal(loc=0, scale = 10)])
#        self.pos = np.array([ np.random.normal(loc=0, scale=self.R*self.kappa*5), np.random.uniform(low= -self.R*kappa, high=self.R*kappa*1000)])
        self.pos = np.array([ np.random.normal(loc=0, scale = self.kappa*5e-4), np.random.uniform(low= -kappa*0.001, high=kappa*0.001)])
        self.F0 = 20e5
        self.l = 1e-9
        self.wavelenght = 5e-3
        self.k  = 2*np.pi/(self.wavelenght)

    def force_surface(self, particles):
        f = np.zeros(2)
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                VanderWaals_force = self.VanderWaals*(2 + 5.32/(100*self.kappa)*distance)/((distance**3)*(1+5.32/(100*self.kappa)*distance))               
                Debye_Huckel_force = self.DLVO*np.exp(-distance)
                hydrophobic_force = self.F0*np.exp(-distance/(self.kappa*self.l))
                f+= (-VanderWaals_force -hydrophobic_force+ Debye_Huckel_force)*dr/distance
        return f
    
    def force_vertical(self, particles):
        lambda_b = 1/(4*np.pi*kb*self.T*epsilon)*e**2
        Yukanawa_coef = lambda_b*self.kappa**2*(np.exp(self.kappa*self.R)/(1+self.kappa*self.R))**2
        f = np.zeros(2)
        for particle in particles:
            if particle.index != self.index:
                dr = particle.pos - self.pos
                distance = np.sqrt(np.dot(dr,dr))
                VanderWaals_force = self.VanderWaals*(2 + 5.32/(100*self.kappa)*distance)/((distance**3)*(1+5.32/(100*self.kappa)*distance))               
                Debye_Huckel_force = 2*Yukanawa_coef*np.exp(-distance)/(distance)
                hydrophobic_force = self.F0*np.exp(-distance/(self.kappa*self.l))
                buyancy = (1.0 - 2.12)*1000*9.81*(4/3*np.pi*(self.R)**3)*np.array([0,1]) #tannis acid density
                f+= (-VanderWaals_force -hydrophobic_force+ Debye_Huckel_force + buyancy)*dr/distance
        return f
    
    def brownian_step_xy(self, force, dt, convection, capillary, v_average):
        f = 1/(2*self.kappa*kb*self.T)*force
        print(f)
        if bool(capillary):
            k = self.k
            A0= 5e-4
            omega = 20
            wave_force = 1/(2*self.kappa*kb*self.T)* self.m*omega**2/(4*np.pi**2)*1/4*(k)**3*A0**2*(v_average*(2*self.D*self.kappa))**2*np.sin(2*k*self.pos[0]/self.kappa)
            print('wf', wave_force)
            print(A0,k)
            print(self.pos/(self.kappa*self.wavelenght))
        else:
            wave_force = 0
        if bool(convection):
            vmax = 300e-6
            convection = np.sin(self.pos[0]*k)*(vmax/(2*self.D*self.kappa))
        
        self.rdot =  f + np.random.normal(scale = np.sqrt(dt), size=2) + convection + wave_force*np.array([1,0])
        print('rdot', self.rdot)
        self.pos += self.rdot*dt
        
    def brownian_step_xz(self, force, dt, convection):
        f = 1/(2*self.kappa*kb*self.T)*force
        rdot =  f + np.random.normal(scale = np.sqrt(dt), size=2) + convection*np.array([np.sin(self.pos[0]*2*np.pi/self.lengh),0])
        self.pos += rdot*dt

        

    
    





        

