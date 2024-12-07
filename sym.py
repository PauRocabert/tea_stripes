import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from euler_brownian import positions
from correlation import analysis
import os 

list_dir = os.listdir('./')
Temp = [310,330,350]
con = [0.01,0.1,1.0]
min_dist = np.zeros(shape=(len(Temp), len(con)))
xy = True
for i,T in enumerate(Temp):
    for j,c in enumerate(con):
        if xy:
            data = f'positions_brownian_xy_c_{c}_T_{T}.npz'
        else:
            data = f'positions_brownian_xz_c_{c}_T_{T}.npz'
        if not(data in list_dir):
            print(data)
            positions(T=T, c=c, xy=xy)
            min_dist[i][j]=analysis(data,c=c, t=T)
        if not(f'animation_c={c}_T_stability={T}.mp4' in list_dir):
            x,y=np.load(data)['x'], np.load(data)['y']
            x_mean, y_mean = np.mean(x[:,-1]), np.mean(y[:,-1])
            x_std, y_std= np.std(x[:,-1]), np.std(y[:,-1])
            x, y = x[(np.abs(x[:,-1]) < x_mean + 2*x_std)*(np.abs(y[:,-1])<y_mean+2*y_std)],y[(np.abs(x[:,-1]) < x_mean + 2*x_std)*(np.abs(y[:,-1])<y_mean+2*y_std)] 
            print(np.shape(x))
            Na = 6.022e23
            I = c*Na*1e3
            radius = 1e-6
            kappa = 2.2912074e-3*np.sqrt(I/T)
            N, t = np.shape(x)
            fig, ax = plt.subplots()
            scater = ax.scatter(x[:,0], y[:,0], s=1)
            ax.set_xlim(np.min(x),np.max(x))
            ax.set_xlabel(r'$x/r$')
            ax.set_ylabel(r'$y/r$')
            ax.grid(visible=True, axis = 'both')
            ax.set_ylim(np.min(y), np.max(y))
            def update(frame):
                scater.set_offsets(np.c_[x[:, frame], y[:, frame]])
                return scater,
            ani = FuncAnimation(fig, update, frames=int(t), blit=True)
            ani.save(f'animation_stationary_c={c}_T={T}.mp4', fps=100) 


        
