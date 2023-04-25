import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def Vol(N, p):
    # number of edges in the complete graph
    E = int(N*(N - 1) / 2)
    # number of nodes in the segmented complete graph
    Nn = N + 2*E
    M_1 = N - 1
    M_2 = 1 / (1 - p)
    return M_1*N + M_2*(Nn - N)

x = np.linspace(0, 1, 102)
x = x[1:-1] # remove end points
vol_10 = np.zeros((100))
vol_5 = np.zeros((100))

for i in range(0, 100):
    vol_10[i] = math.sqrt(Vol(10, x[i]))
    vol_5[i] = math.sqrt(Vol(5, x[i]))

ax = plt.figure().add_subplot(111)
plt.rc('text', usetex=True)
#plt.rcParams['text.usetex'] = True
plt.plot(x, vol_10, label='N = 10')
plt.plot(x, vol_5, label='N = 5')
plt.xlabel('Transition Probability, $p$')
plt.ylabel('$\sqrt{Vol(G_{N,p})}$')
plt.title('$\sqrt{Vol(G_{N,p})}$ vs. Transition Probability')
ax.legend(loc="upper left")
plt.savefig('t_and_vol/vol_vs_p.png')
plt.close()
