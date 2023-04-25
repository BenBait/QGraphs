import matplotlib.pyplot as plt
import numpy as np
plt.rc('text', usetex=True)

# 3 N values for now, 9, 10, 11
# 9 p values 0.1, 0.2, 0.3, 0.4 ...
t_opt = np.zeros((3, 9))

t_opt[0][0] = 9
t_opt[0][1] = 8
t_opt[0][2] = 9
t_opt[0][3] = 9
t_opt[0][4] = 12
t_opt[0][5] = 11
t_opt[0][6] = 14
t_opt[0][7] = 19
t_opt[0][8] = 25
t_opt[1][0] = 10
t_opt[1][1] = 9
t_opt[1][2] = 9
t_opt[1][3] = 9
t_opt[1][4] = 11
t_opt[1][5] = 14
t_opt[1][6] = 14
t_opt[1][7] = 19
t_opt[1][8] = 29
t_opt[2][0] = 10
t_opt[2][1] = 10
t_opt[2][2] = 10
t_opt[2][3] = 11
t_opt[2][4] = 12
t_opt[2][5] = 14
t_opt[2][6] = 18
t_opt[2][7] = 19
t_opt[2][8] = 29

ax = plt.figure().add_subplot(111)
ps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
plt.plot(ps, t_opt[0], label="N=9")
plt.scatter(ps, t_opt[0])
plt.plot(ps, t_opt[1], label="N=10")
plt.scatter(ps, t_opt[1])
plt.plot(ps, t_opt[2], label="N=11")
plt.scatter(ps, t_opt[2])
plt.xlabel('transition probability, $p$')
plt.ylabel('optimal time, $t_{opt}$')
plt.title('Optimal Time vs. Transition Probability')
ax.legend(loc='upper left')
plt.savefig('t_and_vol/t_vs_p.png')
plt.close()
