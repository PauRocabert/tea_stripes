import numpy as np 
import matplotlib.pyplot as plt
density = 2.12e3
R = 1e-4
m = 4/3*np.pi*(R)**3*density 
mu = 8.9e-4
A0 = 1e-7
wavelenght = 5e-3
k = 2*np.pi/(wavelenght)
g = 9.81
omega = 2*np.pi/0.1
a0 = m*g/(mu*omega)*k**2*A0

x0 = wavelenght/7
def sec(x):
    return 1/np.cos(x)
def f(x,t):
    return np.tan(k*x)+sec(k*x)-(np.tan(k*x0)+sec(k*x0))*np.exp(a0*np.cos(omega*t))
def df(x):
    return k*sec(k*x)*(sec(k*x) + np.tan(k*x))

def NR_step(x,t):
    return x-f(x,t)/df(x)


x = np.zeros(shape=1000)
t = np.linspace(0,1,1000)
for n in range(10000):
    x = NR_step(x,t)

plt.plot(t,x)
plt.show()


