import numpy as np
import matplotlib.pyplot as plt

def vz(K=2.682,T=343.15, Ra= 1100.65, p=0, d=1, nu=1):
    k = np.pi*2
    theta = T - 273.15 #temperatura en celsius
    Pr =5000/(theta**2 + 155*theta + 3700)
    Ra = 1100.65
    sigma = 0
    N = 501
    dz = 1/(N-1)
    dz2 = dz**2

    x = np.linspace(0,0.5,N)
    W = np.sin(k*x)
    F = (-k**2 -K**2)*W
    G = (-k**2-K**2-sigma)*F


    epsilon = 10000
    for n in range(1000000):
        if n%1000 ==0:
            print(epsilon)
            print(1/N*np.sqrt(np.dot(W-np.sin(k*x),W-np.sin(k*x))))
        if epsilon < 0.01:
            break
        epsilon_0 = epsilon 
        epsilon =0 
        #condicions de contorn 
        W[0],W[N-1] = 0,0 #W(0)=W(1)=0
        W[1], W[N-2] = 0,0 #W'(0) = W'(1)=0
        G[0], G[N-1] = 0,0
        for i in range(1,N-1):
            if i>1 and i<N-2:
                W_new=1/(2+K**2*dz2)*(W[i+1]+W[i-1]-dz2*F[i])
            else:
                W_new = 0
            F_new=1/(2+(K**2*+sigma)*dz2)*(F[i+1]+F[i-1]-dz2*G[i])
            G_new=1/(2+(K**2*+sigma*Pr)*dz2)*(G[i+1]+G[i-1]+Ra*K**2*dz2*W[i])
            dW = W_new-W[i]
            dF = F_new -F[i]
            dG = G_new - G[i]
            epsilon += np.sqrt(dW**2 + dF**2 + dG**2)
            W[i] = W_new
            F[i] = F_new
            G[i] = G_new
        if epsilon_0 < epsilon:
            break

    plt.plot(x,W)
    np.savez(f'W.npz',W=W)
    plt.show()
vz()


def vel(z):
        vmax = 1000e-6 #m/s
        q0= 7.2569
        q1 = 7.3711
        q2 = 3.664
        vz = vmax*(np.sin(q0*z) - 0.033193*np.cos(q2*z)*np.sinh(q1*z) + 0.033146*np.cosh(q1*z)*np.sin(q2*z))
        return vz