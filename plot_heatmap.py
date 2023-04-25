# data and visualization
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import numpy as np
from scipy.interpolate import griddata
from scipy.signal import argrelextrema
#from scipy.signal import savgol_filter

# scripting capabilities
import os

# build the executable to gather a spectrum of success probabilities
command = "./build.sh 3"
os.system(command)

# number of nodes in the complete graph to try
#num_nodes = [9,10,11]
num_nodes = [5, 10]
# possible transition probabilities, p, to use
ps        = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# possible target nodes, for now just the 1st
targets   = [6]
# keep the maximum time and gamma range consistent
time_max  = 100
#gamma_min = [1, 2, 1, 0.5]
#gamma_max = [3, 4, 2.5, 2]
# FOR N = 10
gamma_step = 0.1
print_in = 0

all_time_series = []
#call the data gathering C executable while iterating through possible params
for n in num_nodes:
    if n == 10:
        gamma_min = [3, 3, 2, 2, 1, 1, 1, 1, 0]
        gamma_max = [8, 8, 6, 6, 4, 4, 4, 3, 3]
    else:
        gamma_min = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        gamma_max = [8, 8, 8, 8, 6, 6, 6, 6, 6]
    i = -1
    for p in ps:
        i+=1 # count for max gamma
        for w in targets:
            params = f"{n} {p} {w} {time_max} {gamma_min[i]} {gamma_max[i]} {gamma_step} {print_in}"
            command = f"./a.out {params}"
            os.system(command)

            dat = np.genfromtxt('success_prob.dat', delimiter=',',skip_header=0)
            x = dat[:,0]
            y = dat[:,1]
            z = dat[:,2]

            # find maximum z value and fix corresponding y value
            max_z_idx = np.argmax(z)
            max_y = y[max_z_idx]

            # filter z values for the fixed y value
            filtered_z = z[y == max_y]
            filtered_x = x[y == max_y]

            #all_time_series.append([filtered_x, filtered_z, max_y])

            # get the optimal time parameter (0th index is first local min)
            '''
            if p <= 0.5:
                t_opt = np.argmax(filtered_z[0:20])
            else:
                t_opt = np.argmax(filtered_z[0:35])
            if p <= 0.93:
                t_opt = np.argmax(filtered_z[0:50])
            else:
                t_opt = np.argmax(filtered_z)
            #t_opt = argrelextrema(filtered_z, np.greater)[0]
            print(f"t_opt = {t_opt}")
            '''

            # plot z as a function of x (pi as a function of t
            plt.plot(filtered_x, filtered_z)
            plt.xlabel('Time, $t$')
            plt.ylabel('Success Probability, $\pi$')
            title_text = f'Success Probability vs. Time for $\gamma$={max_y}, $N$={n}, $p$={p}, $w$={w}'
            plt.title(title_text, fontsize=14)
            #caption = f"N = {n}, p = {p}, w = {w}"
            #plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment='center', fontsize=12)
            plt.tight_layout()
            plt.savefig(f'time_series/map_{n}_{p}_{w}.png')
            plt.close()

'''
            x=np.unique(x)
            y=np.unique(y)
            X,Y = np.meshgrid(x,y)

            Z=z.reshape(len(y),len(x))

            fig, ax = plt.subplots()

            im = ax.pcolormesh(X,Y,Z)
            #im.set_clim(vmin=0, vmax=np.max(Z))  # set colorbar limits
            # Generate evenly spaced tick locations and labels
            vmax = np.max(z)
            yticks = np.linspace(0, vmax, num=5)
            yticklabels = [f'{t:.2f}' for t in yticks]

            # Set the tick locations and labels for the colorbar
            cbar = ax.figure.colorbar(im)
            cbar.set_ticks(yticks)
            cbar.set_ticklabels(yticklabels)
            cbar.set_label('success probability')
            #cbar.ax.set_ylabel('success probability')
            #cbar.ax.set_yticklabels(['0', f'{np.max(Z):.2f}'])
            #fig.colorbar(im)
            ax.set_title('t vs. gamma vs. success probability')
            caption = f"complete graph n = {n}; transition probability = {p}; target node = {w}"
            plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment='center', fontsize=12)
            plt.savefig(f'heatmaps/map_{n}_{p}_{w}.png')
            plt.close()

'''

'''
ax = plt.figure().add_subplot(111)
plt.xlabel('t')
plt.ylabel('success probability')
plt.title('success probability vs. t for gamma = {}'.format(max_y))
caption = f"complete graph n = {n}; target node = {w}"
plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment='center', fontsize=12)

colors = ['r', 'b', 'g', 'y', 'k', 'm', 'c', 'orange', 'purple']
i = 0
for ts in all_time_series:
    filtered_x = ts[0]
    smooth_filtered_z = savgol_filter(ts[1], 51, 3)
    max_y      = ts[2]
    t_opt = argrelextrema(smooth_filtered_z, np.greater)
    print(f"p = {ps[i]}, t_opt_smooth = {t_opt}")

    plt.plot(filtered_x, smooth_filtered_z, colors[i], label=f'p = {ps[i]}, max gamma = {max_y}')
    i+=1

ax.legend(loc="upper right")
plt.savefig(f'time_series/cloud_map_{num_nodes[0]}_{targets[0]}.png')
plt.close()
'''
