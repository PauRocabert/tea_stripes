from particle_class import _particle_
import numpy as np 
import matplotlib.pyplot as plt

def positions(N=100,T=300,c=0.001, xy=True):
    m = 1e-10
    eta = 0.001
    radius = 1e-5
    kb = 1.3806488e-23
    kappa = np.sqrt(c)/0.304*1e9
    tau = 6*np.pi*eta*radius/((kappa)**2*kb*T)
    D = kb*T/(6*np.pi*eta*radius)
    time_steps = 1000
    dt = (10/(tau*time_steps)) #aproximadamanet 10 segon 
    vmax = 3e-4/(2*D*kappa)
    def print_init(particle):
        print('*'*10)
        print('\n')
        print(f'time steps :{time_steps}')
        print(f'tau {tau}')
        print(f'dt: {dt}')
        print(f'concentration: {c} (M)')
        print(f'kappa: {kappa} (m)^-1')
        print(f'particle radius: {radius} (m)')
        print(f'vmax: {vmax}')
        print(f'DLVO_coef: {particle.DLVO}')
        print(f'VanderWaals_coef: {particle.VanderWaals}')

    particles = [_particle_(T, m, c,n) for n in range(N)]
    v_average = [particle.rdot[0] for particle in particles]
    print_init(particles[0])

    x = np.zeros(shape=(N,time_steps))
    y = np.zeros(shape=(N,time_steps))
    for n,particle in enumerate(particles):
        x[n][0], y[n][0] = particle.pos
    plt.scatter(x[:,0]/(particles[0].kappa*radius), y[:,0]/(particles[0].kappa*radius), s=1)
    plt.xlabel(r'$x(\kappa^{-1})$')
    plt.ylabel(r'$y(\kappa^{-1})$')
#    plt.show()
    plt.close()
    if xy:
        for t in range(time_steps):   
            print(t/time_steps)
            forces = [particle.force_surface(particles) for particle in particles]
            for n,particle in enumerate(particles):
                particle.brownian_step_xy(forces[n],dt, 0,0, v_average[n])
                x[n][t], y[n][t] = particle.pos/(kappa*particle.R)
                v_average[n] += 1/(t+1)*particle.rdot[0]
                
        np.savez(f'positions_brownian_xy_c_{c}_T_{T}.npz', x=x, y=y)
    else:
        for t in range(time_steps):
            print(t/time_steps)
            forces = [particle.force_vertical(particles) for particle in particles]
            for n,particle in enumerate(particles):
                particle.brownian_step_xz(forces[n],dt, 0)
                x[n][t], y[n][t] = particle.pos/(kappa*particle.wavelenght)
        np.savez(f'positions_brownian_xz_c_{c}_T_{T}.npz', x=x, y=y)
            

    



        
  
    

        

        
    
        
