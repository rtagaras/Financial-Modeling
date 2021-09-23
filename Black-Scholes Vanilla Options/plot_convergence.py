import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
fig.suptitle("Value as a Function of Number of Trials", fontsize=40)

# Adjusting x-axis ticks
x_steps = 20
x_ticks = [i for i in range(x_steps) if i%4==0]

# Files to load
filenames = ["call_vals", "put_vals", "DC_vals", "DP_vals"]

axis_names = ["Call", "Put", "Digital Call", "Digital Put"]

# Create and configure each subplot, then plot data
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(1, 5):

    # Config
    ax = fig.add_subplot(2, 2, i)
    ax.set_xlabel("Log$_2$(Number of Trials)", fontsize = 40)
    ax.set_ylabel(axis_names[i-1] + " Value", fontsize = 40)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_xticks(x_ticks)

    # Load and plot each datafile
    data = np.loadtxt("./Data/" + filenames[i-1] + ".txt")
    ax.plot(data, marker="o", markersize=10, lw=4)

plt.show()