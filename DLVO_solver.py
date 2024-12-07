import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt 

e = 1.6021766e-19
eta = 1e-3 #Pa
epsilon0 = 8.85541878e-12
epsilon_r = 80.10
epsilon = epsilon_r*epsilon0
Na = 6.022e23
kb = 1.3806488e-23
radius = 1e-5
def f(x,c,T, R):
    I = c*Na*1e3
    phi0 = 0.030 #mV
    kappa = 2.2912074e-3*np.sqrt(I/T)
    x = x*kappa
    z = np.tanh(4*e*phi0/(4*kb*T))
    DLVO = 4*64*kb*T*Na*c*(z**2)/(R*kappa)
    VanderWaals = 0.95e-20*kappa**3/(12*np.pi)
    print(x)
    F0 = 20e5
    l = 1e-9
    VanderWaals_force = VanderWaals/(x**2*(kappa + 5.32/(100)*x))
    Debye_Huckel_force = DLVO*np.exp(-x)
    hydrophobic_force = F0*np.exp(-x/(kappa*l))
    f= -VanderWaals_force +Debye_Huckel_force
    return f
X = np.linspace(1e-7,1e-5,1000)
concentracions = np.array([0.001, 0.01, 0.1,1])
force = np.array([[f(x,c, 300, 1e-5) for x in X] for c in concentracions])
fig, axes = plt.subplots()
for i,con in enumerate(concentracions):
    axes.plot(X, force[i], label=str(-np.log(con)))
axes.legend()
plt.show()
plt.close()




