import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


dat = np.genfromtxt('success_prob.dat', delimiter=',',skip_header=0)
x = dat[:,0]
y = dat[:,1]
z = dat[:,2]

x=np.unique(x)
y=np.unique(y)
X,Y = np.meshgrid(x,y)

Z=z.reshape(len(y),len(x))

fig, ax = plt.subplots()

im = ax.pcolormesh(X,Y,Z)
fig.colorbar(im)
ax.set_title('t vs. gamma vs. success probability')
plt.savefig('map.png')
exit()
#dat = np.genfromtxt('success_prob.dat', delimiter=',',skip_header=0)
# X_dat = dat[:,0]
# Y_dat = dat[:,1]
# Z_dat = dat[:,2]

# # Convert from pandas dataframes to numpy arrays
# X, Y, Z, = np.array([]), np.array([]), np.array([])
# for i in range(len(X_dat)):
        # X = np.append(X, X_dat[i])
        # Y = np.append(Y, Y_dat[i])
        # Z = np.append(Z, Z_dat[i])

# # create x-y points to be used in heatmap
# xi = np.linspace(X.min(), X.max(), 1000)
# yi = np.linspace(Y.min(), Y.max(), 1000)

# # Interpolate for plotting
# zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

# I control the range of my colorbar by removing data
# outside of my range of interest
#zmin = 3
#zmax = 12
#zi[(zi<zmin) | (zi>zmax)] = None

# Create the contour plot
CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow)
                  #vmax=zmax, vmin=zmin)
plt.colorbar()
plt.show()
exit()

#a = np.genfromtxt("success_prob.dat", delimiter=",")
#plt.imshow(a, cmap='hot', interpolation='nearest')
#plt.savefig('map.png')

#pi_interpolated = ((a[:,0], a[:,1]), a[:, 2], method='cubic')o

#CS = plt.contourf(a[:,0], , yi, zi, 15, cmap=plt.cm.rainbow,
#                  vmax=zmax, vmin=zmin)
#plt.colorbar()
#plt.show()
