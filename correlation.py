import numpy as np 
from scipy.stats import linregress
import matplotlib.pyplot as plt

def analysis(data, c, t):
    X,Y=np.load(data)['x'], np.load(data)['y']
    N,T = np.shape(X)
    times = [int(t/100*T) for t in range(100)]
    min_dists =[]
    for time in times:
        x,y = X[:,time],Y[:,time]
        N = np.shape(x)[0]
        distance_matrix = np.zeros(shape=(N,N))
        for n in range(N):
            for m in range(N):
                if n<m:
                    distance_matrix[n][m] = np.sqrt((x[m]-x[n])**2+(y[m]-y[n])**2)
                    distance_matrix[m][n] =  distance_matrix[n][m]
        m = 0
        min_dist = np.average(np.array([np.mean(np.sort(row)[:3]) for row in distance_matrix]))
        min_dists.append(min_dist)
    fig, ax = plt.subplots()
    ax.plot(times, min_dists)
    ax.set_xlabel('time')
    ax.set_ylabel(r'd_min/r')
    plt.savefig(data+'c_vs_time.png', dpi=200)
    plt.close()
"""    for n in range(N):
        ac = np.zeros(shape=len(rs))
        for i in range(len(rs)-1):
            nums =np.sum((distance_matrix[n] > Rs[i])*(distance_matrix[n] < Rs[i+1]))
            ac[i] += nums/(2*rs[i]*dr + dr**2)
        m += linregress(rs, ac).slope
    ac,rs=ac[ac>0], rs[ac>0]
    print(ac)
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(rs, ac)
    plt.savefig(f'g(r,c={c},T={t}).png', dpi=200)
    return min_dist """